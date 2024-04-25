%% format data
rng('shuffle');
seed1 = rng().Seed;
seed2 = seed1;
snr_X = 0; %<<< sweep over 
snr_Y = 0.05; %<<< sweep over
train_time_frac = 0.75; %<<< sweep ovr
subsamp_t = 2;
num_train_inds = 18; %<<< sweep over

noise_alg_X = 'logn'; noise_alg_Y = 'logn'; %<<< fixed
test_length = 40; %<<< fixed
err_tol = 0.2; %<<< fixed
stop_tol = 100; %<<< fixed
toggle_zero_crossing = 1; %<<< fixed

toggle_sim = 1; %<<< doesn't affect alg
toggle_vis = 1; %<<< doesn't affect alg
toggle_view_data = 0;
tol_dd_sim = 10^-10; %<<< doesn't affect alg

% %% set parameters
% phifun_Y = @(t)exp(-5*[1./(1-t.^2)-1]);
% phifun_Y = optTFcos(3,0);
phifun_Y = @(t)(1-t.^2).^9;

tf_Y_params = {'meth','FFT','param',2,'mtmin',3,'subinds',2};%<<< user choice

WENDy_args = {'maxits_wendy',10,...
    'lambdas',10.^linspace(-4,0,50),'alpha',0.02,...
    'ittol',10^-4,'diag_reg',10^-6,'verbose',0};
autowendy = 1; %<<< decision needs to made
tol = 5; %<<< decision needs to made
tol_min = 0.1; %<<< decision needs to made
tol_dd_learn = 10^-8;%<<< decision made

pmax_IC = 4;%<<< decision made 
polys_Y_Yeq = 0:3; %<<< decision needs to made
polys_X_Xeq = 0:2; %<<< decision needs to made
pmax_X_Yeq = 4; %<<< decision made 
pmax_Y_Xeq = 4; %<<< decision made
neg_Y = 0; %<<< decision needs to made
neg_X = 0; %<<< decision needs to made
boolT = @(T)all([min(T,[],2)>=-2 sum(T,2)>=-2 max(T,[],2)<4],2); %<<< decision made

%%% IC map assumed to be polynomial in X, no dependence on Y. Incremental lib for X
%%% Xeq assumed to be power-law in X and Y, fixed lib for X, incremental lib for Y
%%% Yeq assumed to be power-law in X and Y, fixed lib for Y, incremental for X

%% get data
dr = '/home/danielmessenger/Dropbox/Boulder/research/data/dukic collab/';
% load([dr,'FH_feedback.mat']);
warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle')
load([dr,'Gregs_mod_V=0.5.mat'],'Ycell','X','t_epi','custom_tags_X',...
    'yearlength','custom_tags_Y','linregargs_fun_IC','linregargs_fun_Y',...
    'linregargs_fun_X','nstates_X','nstates_Y','W_IC_true','tags_IC_true',...
    'W_Y_true','tags_X_true','tags_Y_true','W_X_true','tags_Ext_X_true','tags_Ext_Y_true',...
    'rhs_IC_true','rhs_Y_true','rhs_X_true','tn','t','Y');

[Y_train,X_train,train_inds,train_time,nstates_X,nstates_Y,X_in,sigma_X,sigma_Y,nX,nY] = ...
    format_data(Ycell,X,t_epi,subsamp_t,train_time_frac,num_train_inds,test_length,snr_X,snr_Y,...
    noise_alg_X,noise_alg_Y,seed1,seed2);

if toggle_view_data==1
    for j=1:nstates_X
        subplot(2,1,j)
        plot(tn,X(:,j),'b-.',t,Y(:,j),'r-o',(train_inds-1)*yearlength,X_train(:,j)*nX(j),'kx','linewidth',3,'markersize',10)
        legend({'X','Y'})
    end
end

% custom_tags_Y = {[1+V zeros(1,nstates_Y-2) 1]}; %<<< decision needs to made
% custom_tags_X = {[-V 0]}; %<<< decision needs to made
linregargs_fun_IC = @(WS){'Aineq',-WS.G{1},'bineq',zeros(size(WS.G{1},1),1)};
% linregargs_fun_Y = @(WS){};
% linregargs_fun_X = @(WS){'Aineq',-WS.G{1},'bineq',zeros(size(WS.G{1},1),1)};
% 
% N = 50; UB = 0.02;
% linregargs_fun_IC = @(WS)enforce_pos_zero(WS);
% linregargs_fun_Y = @(WS)enforce_pos_zero(WS);
% linregargs_fun_X = @(WS)enforce_pos(WS,N,UB);
% 
% linregargs_fun_IC = @(WS){};
linregargs_fun_Y = @(WS){};
linregargs_fun_X = @(WS){};


% sigma_est_X = [];sigma_X; %<<< decision needs to made
% sigma_est_Y = [];sigma_Y; %<<< decision needs to made

