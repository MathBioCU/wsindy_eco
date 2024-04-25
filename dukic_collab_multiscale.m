%% get wsindy_data object
ntraj = length(xcell);
nstates = size(xcell{1},2);
noise_ratio = 0;
rng('shuffle'); rng_seed = rng().Seed; rng(rng_seed);

%%% library
polys = 0:3;
tags = get_tags(polys,[],nstates);
lib = library('tags',tags);

%%% test functions
phifuns = {optTFcos(3,2)};
param_tf = {'meth','param','phifuns'};
param_vals = {{'FFT'},{2},phifuns};

%%% optimizer
optm = WS_opt();
toggle_wendy = 0;

%%% validate
toggle_compare = 1:5;
tol_dd = 10^-12;
test_frac = 1;

%%% get models
Wall = zeros(nstates*length(lib.terms),ntraj);
for n = 1:ntraj 
    Uobj = wsindy_data(xcell{n},t_epi);
    Uobj.addnoise(noise_ratio,'seed',rng_seed);
    
    s = cell(1,length(param_tf)); 
    e = cellfun(@(p)1:length(p),param_vals,'uni',0);
    [s{:}] = ndgrid(e{:});
    s = cell2mat(cellfun(@(ss)ss(:),s,'uni',0));
    tf = [];
    for j=1:size(s,1)
        opt = {Uobj,'subinds',ceil(size(s,1)*1.5)};
        for k=1:size(s,2)
            opt = [opt,{param_tf{k},param_vals{k}{s(j,k)}}];
        end
        tf = [tf;{repelem({testfcn(opt{:})},1,Uobj.nstates)}];
    end
    
    %%% wsindy model
    WS = wsindy_model(repmat(Uobj,size(s,1),1),lib,tf,'coarsen_L',1,'multitest',0);
    
    %%% solve
    if toggle_wendy==0
        [WS,loss_wsindy,its,G,b] = optm.MSTLS(WS,'lambda',10.^linspace(-4,-1,50));
    elseif toggle_wendy==1
        [WS,w_its,res,res_0,CovW] = optm.wendy(WS,'maxits',100,'regmeth','MSTLS');
    elseif toggle_wendy==2
        [WS,loss_wsindy,its,G,b] = optm.MSTLS(WS,'alpha',0.01);
        [WS,w_its,res,res_0,CovW] = optm.wendy(WS,'maxits',20);
    elseif toggle_wendy==3
        [WS,loss_wsindy,lambda,w_its,res,res_0,CovW] = optm.MSTLS_WENDy(WS,'maxits_wendy',20,'lambda',10.^linspace(-3,-1,50));
    end
    
    Wall(:,n) = WS.weights;
end

%%% compare to data
for n=1:length(toggle_compare)
    options_ode_sim = odeset('RelTol',tol_dd,'AbsTol',tol_dd*ones(1,nstates));
    Uobj = wsindy_data(xcell{n},t_epi);
    WS_0 = wsindy_model(Uobj,lib,[]);
    w_test = Wall(:,toggle_compare(n));
    rhs_learned = WS_0.get_rhs('w',w_test);
    x0 = cellfun(@(x)x(1),Uobj.Uobs);
    t_train = Uobj.grid{1}(1:floor(end*test_frac));    
    [t_learned,x_learned]=ode15s(@(t,x)rhs_learned(x),t_train,x0,options_ode_sim);
    figure(1);clf
    for j=1:Uobj.nstates
        subplot(Uobj.nstates,1,j)
        plot(Uobj.grid{1},Uobj.Uobs{j},'b-',t_learned,x_learned(:,j),'--','linewidth',2)
        U_true = Uobj.Uobs{j}(1:min(length(t_learned),floor(end*test_frac)));
        err = norm(x_learned(:,j)-U_true)/norm(U_true);
        title(['rel err=',num2str(err)])
        legend({'data','learned'})
    end
end

%% get interpolation dependence on slow vars
pmax = 10;
tol = 0.01;
p=-1;err = tol*2;
S = mean(Wall~=0,2)>0.5;
tS = find(all(Wall(S,:)~=0));
while and(err>tol,p<pmax)
    p=p+1;
    polys_NZ = 0:p;
    tags_NZ = get_tags(polys_NZ,[],size(NZ,2));
    tags_NZ = [tags_NZ;tags_NZ(2:end,:).*[1 -1];tags_NZ(2:end,:).*[-1 1]];
    A_NZ = cell2mat(arrayfun(@(j)prod(NZ(1:end-1,:).^tags_NZ(j,:),2),1:size(tags_NZ),'uni',0));
    w_NZ = A_NZ(tS,:) \ Wall(S,tS)';
    err = max(vecnorm(A_NZ(tS,:)*w_NZ-Wall(S,tS)')./vecnorm(Wall(S,tS)'));
    if p==pmax
        disp(['tolerance not reached'])
    end
end

W_fcn_cell = cell(length(find(S)),1);
for j=1:length(find(S))
    f = @(X) X(1)*0;
    for k=1:size(tags_NZ,1)
        f = @(X) f(X)+w_NZ(k,j)*prod(X.^tags_NZ(k,:));
    end
    W_fcn_cell{j} = f;
end

WS_0 = wsindy_model(wsindy_data(xcell{1},t_epi),lib,[]);
s = WS.get_supp('w',S);
features = cellfun(@(s)lib.get_fHandles(s),s,'uni',0);
features_W = reshape_cell(W_fcn_cell, cellfun(@(s)length(s),s));
param_map = @(X)cellfun(@(f)cellfun(@(g)g(X),f),features_W,'uni',0);
rhs_xy = @(y,X) rhs_fun(features,param_map(X),y);

%% get initial condition map
IC = cell2mat(cellfun(@(x)x(1,:),xcell,'uni',0));
pmax = 10;
tol = 0.1;
p=-1;err = tol*2;
while and(err>tol,p<pmax)
    p=p+1;
    polys_NZ = 0:p;
    tags_NZ = get_tags(polys_NZ,[],size(NZ,2));
    A_NZ = cell2mat(arrayfun(@(j)prod(NZ(1:end-1,:).^tags_NZ(j,:),2),1:size(tags_NZ),'uni',0));
    IC_NZ = A_NZ \ IC;
    err = max(vecnorm(A_NZ*IC_NZ-IC)./vecnorm(IC));
    if p==pmax
        disp(['tolerance not reached'])
    end
end

IC_map = cell(nstates,1);
for j=1:nstates
    f = @(X) X(1)*0;
    for k=1:size(tags_NZ,1)
        f = @(X) f(X)+IC_NZ(k,j)*prod(X.^tags_NZ(k,:));
    end
    IC_map{j} = f;
end

%% compare micro-model
toggle_compare = [1:5];
for n=1:length(toggle_compare)
    options_ode_sim = odeset('RelTol',tol_dd,'AbsTol',tol_dd*ones(1,nstates));
    rhs_learned = @(y)rhs_xy(y,NZ(toggle_compare(n),:));
    x0 = cellfun(@(f)f(NZ(toggle_compare(n),:)),IC_map);
    t_train = t_epi(1:floor(end*test_frac));
    [t_learned,x_learned]=ode15s(@(t,x)rhs_learned(x),t_train,x0,options_ode_sim);
    figure(1);clf
    for j=1:nstates
        subplot(nstates,1,j)
        plot(t_epi,xcell{toggle_compare(n)}(:,j),'b-',t_learned,x_learned(:,j),'--','linewidth',2)
        U_true = xcell{toggle_compare(n)}(1:min(length(t_learned),floor(end*test_frac)),j);
        err = norm(x_learned(:,j)-U_true)/norm(U_true);
        title(['rel err=',num2str(err)])
        legend({'data','learned'})
    end
end
 
function Cnew = reshape_cell(C,a)
    L = length(C(:));
    if L~=sum(a)
        disp(["Warning: can't reshape"])
    else
        Cnew = cell(length(a),1);
        c = 1;
        for j=1:length(a)
            Cnew{j} = C(c:c+a(j)-1);
            c = c+a(j);
        end
    end
end