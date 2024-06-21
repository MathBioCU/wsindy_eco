dr = '/home/danielmessenger/Dropbox/Boulder/research/data/dukic collab/';
load([dr,'Gregs_mod_V=0.5.mat'],'Ycell','X','t_epi',...
    'yearlength','W_IC_true','tags_IC_true',...
    'W_Y_true','tags_X_true','tags_Y_true','W_X_true','tags_Ext_X_true','tags_Ext_Y_true');

%%% get tags and true sparsity pattern
W = [W_IC_true,W_Y_true,W_X_true];
[w_vec_true,ss] = wendy_param(W);
tags_IC_X=tags_IC_true;
tags_Y_Yeq =tags_Y_true; tags_X_Yeq=tags_X_true;
tags_X_Xeq =tags_Ext_X_true; tags_Y_Xeq=tags_Ext_Y_true;

%%% simulation params
x0 = X(1,:);
num_gen = size(X,1)-1;
sig_tmax = 0;
num_t_epi = 112;
tol_dd_sim = 10^-8;
toggle_zero_crossing=1;
stop_tol = 1000;    

%%% sampled data params
subsamp_t = 2;
train_time_frac = 0.75;
num_train_inds = 16;
test_length = 40;
snr_X = 0.0001;
snr_Y = 0.01;
noise_alg_X = 'logn';
noise_alg_Y = 'logn';
gensamp_seed = 1;
noise_seed = 'shuffle';
[Y_train,X_train,train_inds,train_time,nstates_X,nstates_Y,X_in,sigma_X,sigma_Y,nX,nY] = ...
    format_data(Ycell,X,t_epi,subsamp_t,train_time_frac,num_train_inds,...
    test_length,snr_X,snr_Y,noise_alg_X,noise_alg_Y,gensamp_seed,noise_seed);
X_train = X_train.*nX;

%%% nls params
nls_IO_var = 0.1;
w_vec_IO = w_vec_true.*(1 + sqrt(3)*nls_IO_var*(rand(length(w_vec_true),1)-0.5)*2);
bnds = {[],[]};
options_nls = optimoptions(@lsqnonlin,...
    'Algorithm','Levenberg-Marquardt',...
    'Display','iter',...
    'MaxIterations',1000,'MaxFunctionEvaluations',2000,...
    'OptimalityTolerance',1e-10,'StepTolerance',1e-10);
fun = @(w_vec) obj_fun(w_vec,X_train,train_inds,ss,tags_IC_X,tags_Y_Yeq,tags_X_Yeq,tags_X_Xeq,tags_Y_Xeq,nstates_X,nstates_Y,...
    x0,num_gen,num_t_epi,yearlength,sig_tmax,tol_dd_sim,toggle_zero_crossing,stop_tol);

%%% get NLS sol
tic,
param_nls = lsqnonlin(@(w)fun(w),w_vec_IO,bnds{:},options_nls);
total_time_nls = toc;
err_NLS = norm(w_vec_true-param_nls)/norm(w_vec_true);
disp(err_NLS)

[out,X_pred] = obj_fun(param_nls,X_train,train_inds,ss,tags_IC_X,tags_Y_Yeq,tags_X_Yeq,tags_X_Xeq,tags_Y_Xeq,nstates_X,nstates_Y,...
    x0,num_gen,num_t_epi,yearlength,sig_tmax,tol_dd_sim,toggle_zero_crossing,stop_tol);

n = size(X_pred,1);
err_tol = 0.5;
X_test = X;
cumerr = arrayfun(@(i)norm(vecnorm(X_pred(1:i,:)-X_test(1:i,:),2,2))/norm(vecnorm(X_test(1:i,:),2,2)),(1:n)');
n_err_tol = find(cumerr>err_tol,1);
if isempty(n_err_tol)
    n_err_tol = n-1;
else
    n_err_tol = n_err_tol-1;
end
disp(n_err_tol)

function [out,X_pred] = obj_fun(w_vec,X_train,train_inds,ss,tags_IC_X,tags_Y_Yeq,tags_X_Yeq,tags_X_Xeq,tags_Y_Xeq,nstates_X,nstates_Y,...
    x0,num_gen,num_t_epi,yearlength,sig_tmax,tol_dd_sim,toggle_zero_crossing,stop_tol)
    
    [~,~,~,rhs_IC,rhs_Y,rhs_X] = vec2hybrid_sys(w_vec,ss,...
        tags_IC_X,tags_Y_Yeq,tags_X_Yeq,tags_X_Xeq,tags_Y_Xeq,...
        nstates_X,nstates_Y);

    X_pred = sim_hybrid_fcn(rhs_IC,rhs_Y,rhs_X,x0,nstates_Y,...
                num_gen,num_t_epi,yearlength,sig_tmax,tol_dd_sim,toggle_zero_crossing,stop_tol);

    L = size(X_pred,1);
    if L<max(train_inds)
        out = 1;
    else
        out = 0.5*norm(reshape(X_pred(train_inds,:)-X_train,[],1))^2;
    end

end