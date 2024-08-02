folder = fileparts(which(mfilename)); 
addpath(genpath(folder));
rng('shuffle');

%% data hyperparameters
seed1 = 1;   % seed for random generation selection, can be pre-selected generations, or half-width for peak sampling
seed2 = randi(10^9); % seed for random noise 
% seed2 = rng().Seed; % uncomment to save seed for reproducibility
snr_X = 0.00; % noise level for X
snr_Y = 0.00; % noise level for Y
noise_alg_X = 'logn'; % noise distribution for X
noise_alg_Y = 'AWGN'; % noise distribution for Y

num_train_inds = 25; % number of generations observed / number of gens around each peak (if negative)
train_time_frac = 0.8; % fraction of each generation observed
subsamp_t = 2; % within-generation timescale multiplier

%% algorithmic hyperparameters
toggle_zero_crossing = 1; % halt simulations that are non-positive

phifun_Y = @(t)(1-t.^2).^9; % test function for continuous data
tf_Y_params = {'meth','FFT','param',2,'mtmin',3,'subinds',-3};% test function params

WENDy_args = {'maxits_wendy',5,...
    'lambdas',10.^linspace(-4,0,50),'alpha',0.01,...
    'ittol',10^-4,'diag_reg',10^-6,'verbose',1};
autowendy = 0.95; % confidence level for automatic library incrementation
tol = 5; % default heuristic covariance factor for incrementation, chosen when autowendy = 0.5;
tol_min = 0.1; % lower bound on rel. resid. to increment library, default for covariance severely underestimated
tol_dd_learn = 10^-10; % ODE tolerance for forward solves in computing Y(T)
X_var = [];'true'; % specify variances for discrete vars X in WENDy, [] gives 0, 'true' uses true variances used to generate noise

pmax_IC = 1; % max poly degree for IC solve
polys_Y_Yeq = [1 3]; % Y library for Yeq solve
pmax_X_Yeq = 1; % max poly degree for X terms in Yeq solve
polys_X_Xeq = 0:1; % X library in Xeq solve
pmax_Y_Xeq = 1; % max poly degree for Y terms in Xeq solve
neg_Y = 1; % toggle use negative powers for X terms in Yeq
neg_X = 0; % toggle use negative powers for Y terms in Xeq
boolT = @(T)all([min(T,[],2)>=0 sum(T,2)>=-2 max(T,[],2)<=2 max(T(:,6:7),[],2)==0],2); % restrict poly terms in Yeq
custom_tags_Y = {}; % custom Y tags for Yeq. Example: {[1+V zeros(1,nstates_Y-2) 1]}. Stored with data
custom_tags_X = {}; % custom X tags for Yeq. Stored with data
linregargs_fun_IC = @(WS){}; % addition linear regression args, including constraints, as function of WSINDy model object
linregargs_fun_Y = @(WS){}; % Stored with data
linregargs_fun_X = @(WS){}; % Stored with data
%%% example:
% linregargs_fun_IC = @(WS){'Aineq',-WS.G{1},'bineq',zeros(size(WS.G{1},1),1)}; %%% enforce nonnegative IC map on data

%% post-processing
test_length = 20; % number of generations to test over
err_tol = 0.5; % tol for n_tol= number of generations for which cumulative rel err < tol
stop_tol = 100; % halt simulations if values exceed max observed by this multitude
toggle_sim = 1; % toggle perform diagnostic forward simulation
num_sim = 0; % number of out-of-sample testing simulations
oos_std = 0.5; % std of out-of-sample ICs, uniformly randomly sampled around training IC
toggle_vis = 1; % toggle plot diagnostics
toggle_view_data = 1; % toggle view data before alg runs
tol_dd_sim = 10^-10; % ODE tolerance (abs,rel) for diagnostic sim
yscl = 'linear';

%% format data: Ycell,X,t_epi,
% load('data/forced_FHN.mat','Ycell', 'X', 't_epi', 'yearlength') %%% load in data (minimal vars needed)
% load('~/Desktop/host_multipath.mat') %%% load in data include true model to benchmark

num_gen = size(X,1);
tn = (0:num_gen-1)*yearlength; % discrete time
[Y_train,X_train,train_inds,train_time,nstates_X,nstates_Y,X_in,sigma_X,sigma_Y,nX,nY] = ...
    format_data(Ycell,X,t_epi,subsamp_t,train_time_frac,num_train_inds,test_length,snr_X,snr_Y,...
    noise_alg_X,noise_alg_Y,seed1,seed2);
num_t_epi = length(t_epi{1});
if isequal(X_var,'true')
    X_var = max(sigma_X,0);
end
if toggle_view_data==1 %%% view data
    t = cell2mat(arrayfun(@(i)(i-1)*yearlength+t_epi{i},(1:num_gen-1)','uni',0)); % full continuous time
    Y = cell2mat(Ycell); 
    for j=1:nstates_X
        subplot(nstates_X,1,j)
        plot(tn,X(:,j),'b-.',t,Y(:,j),'r-',(train_inds-1)*yearlength,X_train(:,j)*nX(j),'kx','linewidth',3,'markersize',10)
        legend({'X','Y','I'})
    end
end

%% run alg
tic;
[rhs_IC,W_IC,rhs_Y,W_Y,rhs_X,W_X,lib_Y_IC,lib_X_IC,...
    lib_Y_Yeq,lib_X_Yeq,lib_Y_Xeq,lib_X_Xeq,WS_IC,WS_Yeq,WS_Xeq,...
    CovW_IC,CovW_Y,CovW_X,errs_Yend,loss_IC,loss_Y,loss_X] = ...
    wsindy_eco_fcn(toggle_zero_crossing,stop_tol,phifun_Y,tf_Y_params,WENDy_args,autowendy,...
    tol,tol_min,tol_dd_learn,pmax_IC,polys_Y_Yeq,polys_X_Xeq,pmax_X_Yeq,pmax_Y_Xeq,neg_Y,neg_X,boolT,...
    Y_train,X_train,X_var,train_inds,train_time,t_epi,yearlength,nstates_X,nstates_Y,X_in,nX,nY,...
    custom_tags_X,custom_tags_Y,linregargs_fun_IC,linregargs_fun_Y,linregargs_fun_X);
fprintf('\n runtime: %2.3f \n',toc)

%% compare coefficients
if all(cellfun(@(s)exist(s,'var'),{'W_IC_true','W_Y_true','W_X_true'}))
    coeff_compare;
end

%% sim full system
sim_script;