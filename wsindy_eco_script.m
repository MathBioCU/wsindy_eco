addpath(genpath('wsindy_obj_base'))
rng('shuffle');

%% data hyperparameters
seed1 = 2;   % seed for random generation selection, can just be pre-selected generations
seed2 = 2024;
% seed2 = rng().Seed; % seed for random noise
snr_X = 0.0001; % noise level for X
snr_Y = 0.05; % noise level for Y
noise_alg_X = 'logn'; % noise distribution for X
noise_alg_Y = 'logn'; % noise distribution for Y

num_train_inds = -4; % number of generations observed / number of gens around each peak
train_time_frac = 0.75; % fraction of generations observed
subsamp_t = 2; % within-generation timescale multiplier

%% algorithmic hyperparameters
toggle_zero_crossing = 1; % halt simulations that are non-positive

phifun_Y = @(t)(1-t.^2).^9; % test function for continuous data
tf_Y_params = {'meth','FFT','param',2,'mtmin',3,'subinds',-3};% test function params

WENDy_args = {'maxits_wendy',5,...
    'lambdas',10.^linspace(-4,-1,50),'alpha',0.01,...
    'ittol',10^-4,'diag_reg',10^-6,'verbose',1};
autowendy = 0.95; % increment library approximate confidence interval with this confidence level 
tol = 5; % default heuristic increment, chosen when autowendy = 0.5;
tol_min = 0.1; % lower bound on rel. residual to increment library, in case covariance severely underestimated
tol_dd_learn = 10^-10;% ODE tolerance for forward solves in computing Y(T)
X_var = 'true';

pmax_IC = 4; % max poly degree for IC solve
polys_Y_Yeq = 0:3; % Y library for Yeq solve
pmax_X_Yeq = 4; % max poly degree for X terms in Yeq solve
polys_X_Xeq = 0:2; % X library in Xeq solve
pmax_Y_Xeq = 3; % max poly degree for Y terms in Xeq solve
neg_Y = 0; % toggle use negative powers for X terms in Yeq
neg_X = 0; % toggle use negative powers for Y terms in Xeq
boolT = @(T)all([min(T,[],2)>=-2 sum(T,2)>=-2 max(T,[],2)<4],2); % restrict poly terms in Yeq
custom_tags_Y = {}; % custom Y tags for Yeq. Example: {[1+V zeros(1,nstates_Y-2) 1]}
custom_tags_X = {}; % custom X tags for Yeq
linregargs_fun_IC = @(WS){}; % addition linear regression args, including constraints, as function of WSINDy model object
linregargs_fun_Y = @(WS){};
linregargs_fun_X = @(WS){};
%%% example:
% linregargs_fun_IC = @(WS){'Aineq',-WS.G{1},'bineq',zeros(size(WS.G{1},1),1)}; %%% enforce nonnegative IC map on data

%% post-processing
test_length = 40; % number of generations to test over
err_tol = 0.5; % tol for n_tol= number of generations for which cumulative rel err < tol
stop_tol = 100; % halt simulations if values exceed max observed by this multitude
toggle_sim = 1; % toggle perform diagnostic forward simulation
num_sim = 0; % number of out-of-sample testing simulations
oos_std = 0.2; % std of out-of-sample ICs, uniformly randomly sampled around training IC
toggle_vis = 1; % toggle plot diagnostics
toggle_view_data = 0; % toggle view data before alg runs
tol_dd_sim = 10^-10; % ODE tolerance (abs,rel) for diagnostic sim

%% get data
dr = '/home/danielmessenger/Dropbox/Boulder/research/data/dukic collab/';
% load([dr,'FH_feedback.mat']);
warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle')
load([dr,'Gregs_mod_V=0.5.mat'],'Ycell','X','t_epi','custom_tags_X',...
    'yearlength','custom_tags_Y','linregargs_fun_IC','linregargs_fun_Y',...
    'linregargs_fun_X','nstates_X','nstates_Y','W_IC_true','tags_IC_true',...
    'W_Y_true','tags_X_true','tags_Y_true','W_X_true','tags_Ext_X_true','tags_Ext_Y_true',...
    'rhs_IC_true','rhs_Y_true','rhs_X_true','tn','t','Y','sig_tmax');
% gregs_evoMod;

%% format data
[Y_train,X_train,train_inds,train_time,nstates_X,nstates_Y,X_in,sigma_X,sigma_Y,nX,nY] = ...
    format_data(Ycell,X,t_epi,subsamp_t,train_time_frac,num_train_inds,test_length,snr_X,snr_Y,...
    noise_alg_X,noise_alg_Y,seed1,seed2);
num_gen = size(X,1)-1;
num_t_epi = length(t_epi{1});
if isequal(X_var,'true')
    X_var = max(sigma_X,0);
end
if toggle_view_data==1 %%% view data
    for j=1:nstates_X
        subplot(2,1,j)
        plot(tn,X(:,j),'b-.',t,Y(:,j),'r-',(train_inds-1)*yearlength,X_train(:,j)*nX(j),'kx','linewidth',3,'markersize',10)
        legend({'X','Y','I'})
    end
    pause
end

%% run alg
tic;
[rhs_IC,W_IC,rhs_Y,W_Y,rhs_X,W_X,lib_Y_IC,lib_X_IC,...
    lib_Y_Yeq,lib_X_Yeq,lib_Y_Xeq,lib_X_Xeq,WS_IC,WS_Yeq,WS_Xeq,...
    CovW_IC,CovW_Y,CovW_X,errs_Yend] = ...
    wsindy_eco_fcn(toggle_zero_crossing,stop_tol,phifun_Y,tf_Y_params,WENDy_args,autowendy,...
    tol,tol_min,tol_dd_learn,pmax_IC,polys_Y_Yeq,polys_X_Xeq,pmax_X_Yeq,pmax_Y_Xeq,neg_Y,neg_X,boolT,...
    Y_train,X_train,X_var,train_inds,train_time,t_epi,yearlength,nstates_X,nstates_Y,X_in,nX,nY,...
    custom_tags_X,custom_tags_Y,linregargs_fun_IC,linregargs_fun_Y,linregargs_fun_X);
fprintf('\n runtime: %2.3f \n',toc)

%% compare coefficients
coeff_compare;

%% sim full system
sim_script;

%% extra 
% phifun_Y = @(t)exp(-5*[1./(1-t.^2)-1]);
% phifun_Y = optTFcos(3,0);
% phifun_Y = 'delta'; % test function for continuous data
% tf_Y_params = {'meth','direct','param',1,'mtmin',1};% test function params
