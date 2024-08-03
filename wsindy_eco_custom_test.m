folder = fileparts(which(mfilename)); 
addpath(genpath(folder));
rng('shuffle');

%% data hyperparameters
seed1 = 2;   % seed for random generation selection, can be pre-selected generations, or half-width for peak sampling
seed2 = 1;randi(10^9); % seed for random noise 
% seed2 = rng().Seed; % uncomment to save seed for reproducibility
snr_X = 0.00; % noise level for X
snr_Y = 0.00; % noise level for Y
noise_alg_X = 'logn'; % noise distribution for X
noise_alg_Y = 'logn'; % noise distribution for Y

num_train_inds = -5; % number of generations observed / number of gens around each peak (if negative)
train_time_frac = 0.75; % fraction of each generation observed
subsamp_t = 3; % within-generation timescale multiplier
toggle_scale = 1;

%% algorithmic hyperparameters
toggle_zero_crossing = 0; % halt simulations that are non-positive

phifun_Y = @(t)(1-t.^2).^9; % test function for continuous data
tf_Y_params = {'meth','direct','param',3,'mtmin',3,'subinds',-3};% test function params

WENDy_args = {'maxits_wendy',5,...
    'lambdas',10.^linspace(-3,-1,40),'alpha',0.01,...
    'ittol',10^-4,'diag_reg',10^-6,'verbose',1};
autowendy = 0.95; % confidence level for automatic library incrementation
tol = 5; % default heuristic covariance factor for incrementation, chosen when autowendy = 0.5;
tol_min = 0.1; % lower bound on rel. resid. to increment library, default for covariance severely underestimated
tol_dd_learn = 10^-10; % ODE tolerance for forward solves in computing Y(T)
X_var = [];'true'; % specify variances for discrete vars X in WENDy, [] gives 0, 'true' uses true variances used to generate noise

pmax_IC = 3; % max poly degree for IC solve
polys_Y_Yeq = 0:3; % Y library for Yeq solve
pmax_X_Yeq = 3; % max poly degree for X terms in Yeq solve
polys_X_Xeq = 0:2; % X library in Xeq solve
pmax_Y_Xeq = 3; % max poly degree for Y terms in Xeq solve
neg_Y = 0; % toggle use negative powers for X terms in Yeq
neg_X = 0; % toggle use negative powers for Y terms in Xeq
boolT = @(T)all([max(T,[],2)<=2 max(T(:,6:end),[],2)==0],2); % restrict poly terms in Yeq
% boolTL = [repmat({@(T,L) min(T(:,min(find(L.ftag),end)),[],2)>0},1,5),{[]}];
boolTL = repmat({@(T,L) min(T(:,min(find(L.ftag),end)),[],2)>0},1,6);
custom_tags_Y = {[],[],[],[],[],[1 1 0 1 0 0]}; % custom Y tags for Yeq. Example: {[1+V zeros(1,nstates_Y-2) 1]}. Stored with data
custom_tags_X = {}; % custom X tags for Yeq. Stored with data
linregargs_fun_IC = @(WS){}; % addition linear regression args, including constraints, as function of WSINDy model object
linregargs_fun_Y = @(WS){}; % Stored with data
linregargs_fun_X = @(WS){}; % Stored with data
%%% example:
% linregargs_fun_X = @(WS){'Aineq',-WS.G{1},'bineq',zeros(size(WS.G{1},1),1)}; %%% enforce nonnegative IC map on data

%% post-processing
test_length = 20; % number of generations to test over
err_tol = 0.5; % tol for n_tol= number of generations for which cumulative rel err < tol
stop_tol = 100; % halt simulations if values exceed max observed by this multitude
toggle_sim = 1; % toggle perform diagnostic forward simulation
num_sim = 5; % number of out-of-sample testing simulations
oos_std = 0.5; % std of out-of-sample ICs, uniformly randomly sampled around training IC
toggle_vis = 1; % toggle plot diagnostics
toggle_view_data = 1; % toggle view data before alg runs
tol_dd_sim = 10^-12; % ODE tolerance (abs,rel) for diagnostic sim
yscl = 'log';

%% format data: Ycell,X,t_epi,
% load('data/forced_FHN.mat','Ycell', 'X', 't_epi', 'yearlength') %%% load in data (minimal vars needed)
% load('~/Desktop/host_multipath.mat') %%% load in data include true model to benchmark
% load('~/Desktop/host_multipath_davies.mat') %%% load in data include true model to benchmark
load('~/Desktop/host_multipath_6-3_c.mat','Ycell','X','t_epi','yearlength',...
    'W_IC_true','W_X_true','W_Y_true','tags_Ext_X_true','tags_Ext_Y_true','tags_IC_true','tags_X_true','tags_Y_true',...
    'rhs_Y_true','rhs_X_true','rhs_IC_true') %%% load in data include true model to benchmark

% load('~/Desktop/host_multipath_6-3_d_steady_state_1path.mat','Ycell','X','t_epi','yearlength',...
%     'W_IC_true','W_X_true','W_Y_true','tags_Ext_X_true','tags_Ext_Y_true','tags_IC_true','tags_X_true','tags_Y_true',...
%     'rhs_Y_true','rhs_X_true','rhs_IC_true') %%% load in data include true model to benchmark

num_gen = 50;

num_gen = min(num_gen,size(X,1));
X = X(1:num_gen,:);
Ycell = Ycell(1:num_gen-1);
t_epi = t_epi(1:num_gen-1);

tn = (0:num_gen-1)*yearlength; % discrete time
[Y_train,X_train,train_inds,train_time,nstates_X,nstates_Y,X_in,sigma_X,sigma_Y,nX,nY] = ...
    format_data(Ycell,X,t_epi,subsamp_t,train_time_frac,num_train_inds,test_length,snr_X,snr_Y,...
    noise_alg_X,noise_alg_Y,seed1,seed2,toggle_scale);
num_t_epi = length(t_epi{1});
if isequal(X_var,'true')
    X_var = max(sigma_X,0);
end
if toggle_view_data==1 %%% view data
    t = cell2mat(arrayfun(@(i)(i-1)*yearlength+t_epi{i},(1:num_gen-1)','uni',0)); % full continuous time
    Y = cell2mat(Ycell); 
    for j=1:nstates_Y
        subplot(nstates_Y,1,j)
        try 
            plot(tn,X(:,j),'b-.',t,Y(:,j),'r-',(train_inds-1)*yearlength,X_train(:,j)*nX(j),'kx','linewidth',3,'markersize',10)
            legend({'X','Y','I'})
        catch
            plot(t,Y(:,j),'r-','linewidth',3,'markersize',10)
            legend({'Y'})
        end
        set(gca,'Yscale',yscl)
        grid on
    end
end

%% run alg

% [rhs_IC,W_IC,rhs_Y,W_Y,rhs_X,W_X,...
%     lib_Y_IC,lib_X_IC,lib_Y_Yeq,lib_X_Yeq,lib_Y_Xeq,lib_X_Xeq,...
%     WS_IC,WS_Yeq,WS_Xeq,...
%     CovW_IC,CovW_Y,CovW_X,...
%     Y_ns,errs_Yend,...
%     loss_IC,loss_Y,loss_X]= ...
%     wsindy_eco_fcn(toggle_zero_crossing,stop_tol,phifun_Y,tf_Y_params,WENDy_args,autowendy,tol,tol_min,...
%         tol_dd_learn,pmax_IC,polys_Y_Yeq,polys_X_Xeq,pmax_X_Yeq,pmax_Y_Xeq,neg_Y,neg_X,boolT,boolTL,...
%         Y_train,X_train,X_var,train_inds,train_time,t_epi,yearlength,nstates_X,nstates_Y,X_in,nX,nY,...
%         custom_tags_X,custom_tags_Y,linregargs_fun_IC,linregargs_fun_Y,linregargs_fun_X);
% 

tic;

if isempty(X_var)
    X_var = X_train*0;
end

addpath(genpath('wsindy_obj_base'))
%%% get wsindy_data 
Uobj_Y = cellfun(@(x,t)wsindy_data(x,t),Y_train,train_time);
Uobj_tot = arrayfun(@(i)...
    wsindy_data([[Uobj_Y(i).Uobs{:}] repmat(X_train(X_in(i),:),Uobj_Y(i).dims,1)],train_time{i}),(1:length(Uobj_Y))');
foo = arrayfun(@(U)U.estimate_sigma('set',true),Uobj_tot,'uni',0);
for j=1:length(Uobj_tot)
    Uobj_tot(j).sigmas(nstates_Y+1:end) = num2cell(X_var(X_in(j),:));
end

IC = zeros(max(train_inds),nstates_Y+nstates_X);
IC(train_inds,nstates_Y+1:end) = X_train;
IC(train_inds(X_in),1:nstates_Y) = cell2mat(arrayfun(@(U)cellfun(@(x)x(1,:),U.Uobs),Uobj_Y,'uni',0));
Uobj_IC = wsindy_data(IC,0:max(train_inds)-1);
if autowendy>0
    S = arrayfun(@(U)cell2mat(U.estimate_sigma).^2,Uobj_Y,'uni',0);
    IC_cov = IC*0; 
    IC_cov(train_inds(X_in),1:nstates_Y) = cell2mat(S);
    IC_cov(train_inds(X_in),nstates_Y+1:end) = X_var(X_in,:).^2;
    Uobj_IC.R0 = spdiags(IC_cov(:),0,numel(IC_cov),numel(IC_cov));
end
E = eye(nstates_Y+nstates_X);

%% get IC_map
lib_Y_IC = library('tags',zeros(1,nstates_Y));
lib_X_IC = library();
tf_IC = testfcn(Uobj_IC,'meth','direct','param',0,'phifuns','delta','mtmin',0,'subinds',train_inds(X_in));
lhs_IC = arrayfun(@(i)term('ftag',E(i,:)),(1:nstates_Y)','uni',0);
[rhs_IC,W_IC,WS_IC,lib_X_IC,loss_IC,lambda_IC,W_its_IC,res_IC,res_0_IC,CovW_IC] = ...
    hybrid_MI(pmax_IC,lib_Y_IC,lib_X_IC,nstates_Y,nstates_X,Uobj_IC,tf_IC,lhs_IC,...
                WENDy_args,linregargs_fun_IC,autowendy,tol,tol_min,nY,nX);
rhs_IC = @(X) rhs_IC(zeros(nstates_Y,1),X(:));

%% get parametric small scale model
tf_Yeq = arrayfun(@(U)...
    arrayfun(@(i)testfcn(U,tf_Y_params{:},'phifuns',phifun_Y,'stateind',i),1:U.nstates,'uni',0),...
    Uobj_Y,'uni',0);
lhs_Yeq = arrayfun(@(i)term('ftag',E(i,:),'linOp',1),(1:nstates_Y)','uni',0);
lib_Y_Yeq = arrayfun(@(i)library('nstates',nstates_Y),1:nstates_Y);

if ~isempty(custom_tags_Y)
    if length(custom_tags_Y)==1
        custom_tags_Y = repmat(custom_tags_Y,1,nstates_Y);
    end
    arrayfun(@(L,i)L.add_terms(custom_tags_Y{i}),lib_Y_Yeq,1:nstates_Y);
end
for i=1:nstates_Y
    lib_Y_Yeq(i).add_tags('polys',polys_Y_Yeq,'neg',neg_Y,'boolT',boolT,'boolTL',boolTL{i},'lhs',lhs_Yeq{i}); % library
end
lib_X_Yeq = library('tags',custom_tags_X);
[rhs_Y,W_Y,WS_Yeq,lib_X_Yeq,loss_Y,lambda_Y,W_its_Y,res_Y,res_0_Y,CovW_Y] = ...
    hybrid_MI(pmax_X_Yeq,lib_Y_Yeq,lib_X_Yeq,nstates_Y,nstates_X,Uobj_tot,tf_Yeq,lhs_Yeq,...
                WENDy_args,linregargs_fun_Y,autowendy,tol,tol_min,nY,nX);

%% get Y(T)
X_sub = find(diff(train_inds)==1);
subinds = train_inds(X_sub);
Yend = zeros(length(X_sub),nstates_Y);
X_n = X_train(X_sub,:);
errs_Yend = []; Y_new = {}; Y_ns = cell(length(X_sub),2);
for n=1:length(X_sub)
    options_ode_sim = odeset('RelTol',tol_dd_learn,'AbsTol',tol_dd_learn*ones(1,nstates_Y),'Events',@(T,Y)myEvent(T,Y,stop_tol*max(max(cell2mat(Y_train))),toggle_zero_crossing));
    rhs_learned = @(y)rhs_Y(y,X_n(n,:).*nX);
    Y0 = rhs_IC(X_n(n,:).*nX);
    t_train = t_epi{train_inds(X_sub(n))}([1 end]);
    [t_n,Y_n,TE]=ode15s(@(t,x)rhs_learned(x),t_train,Y0,options_ode_sim);
    Y_ns(n,:) = {Y_n,t_n};
    if ismember(X_sub(n),X_in)
        try
            Y_new = [Y_new;{interp1(t_n,Y_n,train_time{X_in==X_sub(n)})./nY}];
            errs_Yend = [errs_Yend;std(Y_new{end}-Y_train{X_in==X_sub(n)})];
            % figure(1)
            % for i=1:nstates_Y
            %     subplot(nstates_Y,1,i)
            %     plot(t_n,Y_n(:,i)/nY(i),train_time{X_in==X_sub(n)},Y_train{X_in==X_sub(n)}(:,i)); 
            %     legend({'learn','data'})
            % end
            % drawnow
        catch
            Y_new = [Y_new;{NaN}];
            errs_Yend = [errs_Yend;max(abs(Y_n),'omitnan')];
        end
    % figure(1)
    % for i=1:nstates_Y
    %     subplot(nstates_Y,1,i)
    %     plot(t_n,Y_n(:,i)/nY(i),train_time{X_in==X_sub(n)},Y_train{X_in==X_sub(n)}(:,i),t_epi{train_inds(X_sub(n))},Ycell{train_inds(X_sub(n))}(:,i)/nY(i)); 
    %     legend({'learn','data','true'})
    % end
    % drawnow
    end
    Yend(n,:) = Y_n(end,:);
end
Yend = Yend./nY;

%% get large scale model
X_Yend = zeros(max(train_inds),nstates_X+nstates_Y);
X_Yend(train_inds,1:nstates_X) = X_train;
% Tends = cellfun(@(t)t(end),t_epi(subinds));
% Yend = cell2mat(arrayfun(@(i)interp1(Y_ns{i,2},Y_ns{i,1},Tends(i))./nY,(1:length(subinds))','uni',0));
X_Yend(subinds,nstates_X+1:end) = Yend;
Uobj_X_Yend = wsindy_data(X_Yend,(0:max(train_inds)-1)*yearlength);
if autowendy>0
    Xn_cov = X_Yend*0; 
    errs_Yend = fillmissing(interp1(train_inds(X_in),errs_Yend,subinds,'linear'),'linear');
    Xn_cov(subinds,1:nstates_X) = X_var(X_sub,:).^2;
    Xn_cov(subinds,nstates_X+1:end) = errs_Yend.^2;
    Uobj_X_Yend.R0 = spdiags(Xn_cov(:),0,numel(Xn_cov),numel(Xn_cov));
end
lib_X_Xeq = library('tags',get_tags(polys_X_Xeq,[],nstates_X));
lib_Y_Xeq = library();
tf_X = testfcn(Uobj_X_Yend,'meth','direct','param',1/2,'phifuns','delta','mtmin',0,'subinds',subinds);
lhs = arrayfun(@(i)term('ftag',E(i,:),'linOp',1),(1:nstates_X)','uni',0);
[rhs_X,W_X,WS_Xeq,lib_Y_Xeq,loss_X,lambda_X,w_its,res_X,res_0_X,CovW_X] = ...
    hybrid_MI(pmax_Y_Xeq,lib_X_Xeq,lib_Y_Xeq,nstates_X,nstates_Y,...
    Uobj_X_Yend,tf_X,lhs,WENDy_args,linregargs_fun_X,autowendy,tol,tol_min,nX,nY);

fprintf('\n runtime: %2.3f \n',toc)

%% compare coefficients
if all(cellfun(@(s)exist('s','var'),{'W_IC_true','W_Y_true','W_X_true'}))
    coeff_compare;
end

%% sim full system
sim_script;