%%% get random coefficients
n = 1;
W_IC_param = mvnrnd(cell2mat(WS_IC.get_params),CovW_IC,n);
W_Y_param = mvnrnd(cell2mat(WS_Yeq.get_params),CovW_Y,n);
W_X_param = mvnrnd(cell2mat(WS_Xeq.get_params),CovW_X,n);

W_temp = W_X;
ind = 0;
x = W_X{1};
s = find(x);
x(s) = W_X_param(ind+(1:length(s)));
ind = ind + length(s);
W_temp{1} = x;

%%

[rhs_IC_true,W_IC,tags_X_IC] = shorttime_map(W_IC_true,library('tags',zeros(1,nstates_Y)),tags_IC_true,ones(1,nstates_Y),ones(1,nstates_X));
rhs_IC_true = @(X) rhs_IC_true(zeros(nstates_Y,1),X(:));
[rhs_Y_true,W_Y,tags_Y_Yeq] = shorttime_map(W_Y_true,tags_Y_true,tags_X_true,ones(1,nstates_Y),ones(1,nstates_X));
[rhs_X_true,W_X] = longtime_map(W_X_true,tags_Ext_X_true,tags_Ext_Y_true,ones(1,nstates_X),ones(1,nstates_Y));


[rhs,W] = shorttime_map(W,lib_state,lib_param,scale_state,scale_param);
[rhs,W] = shorttime_map(W,lib_state,lib_param,scale_state,scale_param);
[rhs,W] = shorttime_map(W,lib_state,lib_param,scale_state,scale_param);

[X_pred,Ycell_pred,Y_pred,t_pred,tn_pred,t_epi_pred] = sim_hybrid_fcn(...
    rhs_IC,rhs_Y,rhs_X,x0,nstates_Y,...
    num_gen,num_t_epi,yearlength,sig_tmax,...
    tol_dd_sim,toggle_zero_crossing,...
    stop_tol*max(max(cell2mat(Y_train)./nY)));
