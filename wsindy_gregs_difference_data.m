%% get wsindy_data object
nx = [1 1]./mean(cell2mat(cellfun(@(x) mean(x(:)),xcell,'uni',0)));
xred = cellfun(@(x)x.*nx,xcell,'uni',0);
ntraj = length(xred);
tred = t;
Uobj = arrayfun(@(i)wsindy_data(xred{i},tred(:)),(1:ntraj)');
Uobj.trimstart(500);

nstates = Uobj.nstates;
M = Uobj.dims;

noise_ratio = 0;
% rng('shuffle')
rng_seed = rng().Seed; rng(rng_seed);
Uobj.addnoise(noise_ratio,'seed',rng_seed);

%% get lib tags

polys = 0:3;
trigs = [];
tags_1 = get_tags(polys,[],nstates);
tags_1 = unique([tags_1;tags_1.*[-1 1];tags_1.*[1 -1];tags_1.*[-1 -1]],'rows');
tags_1 = tags_1(min(tags_1,[],2)>-3,:);
tags_2 = get_tags([],trigs,nstates);
tags = prodtags(tags_1,tags_2);
lib = library('tags',tags);
tags_comp = tags_1(and(max(tags_1,[],2)<2,min(tags_1,[],2)>0),:);
tags_comp = mat2cell(tags_comp,ones(size(tags_comp,1),1),size(tags_comp,2));
lib.complib({@(x)x./(0.001+x.^2)},tags_comp);
lib.complib({@(x)1./(0.001+x.^2)},{[0 1];[1 0];[1 1]});

%% get test function

phifuns = {optTFcos(3,0),optTFcos(3,1)};
param_tf = {'meth','param','phifuns'};
% param_vals = {{'direct'},num2cell([9 12 15 21]),phifuns};
param_vals = {{'FFT'},num2cell([2]),phifuns};

s = cell(1,length(param_tf));
e = cellfun(@(p)1:length(p),param_vals,'uni',0);
[s{:}] = ndgrid(e{:});
s = cell2mat(cellfun(@(ss)ss(:),s,'uni',0));
tf = [];
for j=1:size(s,1)
    for ll=1:ntraj
        opt = {Uobj(ll),'subinds',ceil(size(s,1)*1.5)};
        for k=1:size(s,2)
            opt = [opt,{param_tf{k},param_vals{k}{s(j,k)}}];
        end
        tf = [tf;{repelem({testfcn(opt{:})},1,nstates)}];
    end
end

%% build WSINDy linear system

WS = wsindy_model(repmat(Uobj,size(s,1),1),lib,tf,'coarsen_L',1,'multitest',0);

%% solve

optm = WS_opt();
toggle_wendy = 0;

if toggle_wendy==0
    [WS,loss_wsindy,its,G,b] = optm.MSTLS(WS,'lambda',10.^linspace(-3,-1,50));
elseif toggle_wendy==1
    [WS,w_its,res,res_0,CovW] = optm.wendy(WS,'maxits',100,'regmeth','MSTLS');
elseif toggle_wendy==2
    [WS,loss_wsindy,its,G,b] = optm.MSTLS(WS,'alpha',0.01);
    [WS,w_its,res,res_0,CovW] = optm.wendy(WS,'maxits',20);
elseif toggle_wendy==3
    [WS,loss_wsindy,lambda,w_its,res,res_0,CovW] = optm.MSTLS_WENDy(WS,'maxits_wendy',20,'lambda',10.^linspace(-3,-1,50));
end

%% simulate learned and true reduced systems

toggle_compare = 1:ntraj;
if ~isempty(toggle_compare)
    w_plot = WS.weights;
    rhs_learned = WS.get_rhs;
    tol_dd = 10^-12;
    options_ode_sim = odeset('RelTol',tol_dd,'AbsTol',tol_dd*ones(1,nstates));
    for i=toggle_compare

        tol_dd = 10^-10;
        x0 = cellfun(@(x)x(1),Uobj(i).Uobs);
        tcap =  50;
        t_data = Uobj(i).grid{1}(1:tcap);
        t_train = linspace(t_data(1),t_data(end),1000);
        x_data = cell2mat(cellfun(@(U)U(1:tcap),Uobj(i).Uobs,'uni',0));%xred{i}(1:tcap,:);
        minfun = @(s) scale_fun(WS,s,tol_dd,x0,t_train,t_data,x_data);
        optimoptions(@fmincon,'MaxFunctionEvaluations',100);
        % s = fmincon(@(s) scale_fun(WS,s,tol_dd,x0,t_train,t_data,x_data),0.9,[],[],[],[],0.8,1);
        s = 1;%0.75;
        rhs_learned = WS.get_rhs('w',WS.weights*s);

        t_train = linspace(Uobj(i).grid{1}(1),Uobj(i).grid{1}(end),2000);    
        [t_learned,x_learned]=ode15s(@(t,x)rhs_learned(x),t_train,x0,options_ode_sim);
        figure(1);clf
        for j=1:nstates
            subplot(nstates,1,j)
            plot(Uobj(i).grid{1},Uobj(i).Uobs{j},'b-',t_learned,x_learned(:,j),'--','linewidth',2)
            try
                title(['rel err=',num2str(norm(x_learned(:,j)-Uobj(i).Uobs{j})/norm(Uobj(i).Uobs{j}))])
            catch
            end
            legend({'data','learned'})
        end
    end
end

function out = scale_fun(WS,s,tol_dd,x0,t_train,t_data,x_data)
    rhs_learned = WS.get_rhs('w',WS.weights*s);
    options_ode_sim = odeset('RelTol',tol_dd,'AbsTol',tol_dd*ones(1,length(x0)));
    [t_train,x_learned]=ode15s(@(t,x)rhs_learned(x),t_train,x0,options_ode_sim);
    x_learned = interp1(t_train,x_learned,t_data);
    out = norm(x_data(:)-x_learned(:));
    % plot(t_data,x_data,t_data,x_learned,'--')
    % drawnow
end