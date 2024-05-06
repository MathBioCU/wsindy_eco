%% red-text model terms

%%% load data
wsindy_eco_inputs; tic;

%%% specify model
V = 0.5;
tags_IC_X = [[1 0];[0 1];[0 0]]; %%% model libraries
tags_Y_Yeq = [[1+V 1];[1 1];[2 1];[0 1]];
tags_X_Yeq = [[-V 0];[0 0]];
tags_X_Xeq = [[1 0];[0 1];[0 0]];
tags_Y_Xeq = [[1 0];[0 0]];
W_IC_sparse = {[1 -4.93e-7 -1.03e-7],[5.31e-3 1.01 -7.08e-5]}; %%% sparsity patterns
W_Y_sparse = {[[0.812 -1.94];[0 0.099];[3.21 0];[0 0]],...
    [[-1.6 0];[0.142 0];[12.5 0];[0 -0.011]]};
W_X_sparse = {[[0.021 0];[0 0];[5.49 0]],[[0 6.97];[0 0.3];[-6.98 0]]};

%% extract model terms from wsindy_eco run

WENDy_args = {'maxits',10,'ittol',10^-4,'diag_reg',10^-6,'verbose',1};
[rhs_IC,W_IC,rhs_Y,W_Y,rhs_X,W_X,...
    lib_Y_IC,lib_X_IC,lib_Y_Yeq,lib_X_Yeq,lib_Y_Xeq,lib_X_Xeq,...
    WS_IC,WS_Yeq,WS_Xeq,...
    CovW_IC,CovW_Y,CovW_X,...
    errs_Yend]= ...
    wendy_eco_fcn(toggle_zero_crossing,stop_tol,phifun_Y,tf_Y_params,WENDy_args,autowendy,tol_dd_learn,...
        Y_train,X_train,X_var,train_inds,train_time,t_epi,yearlength,nstates_X,nstates_Y,X_in,nX,nY,...
        tags_IC_X,tags_Y_Yeq,tags_X_Yeq,tags_X_Xeq,tags_Y_Xeq,W_IC_sparse,W_Y_sparse,W_X_sparse);

%% display models

IC_mod = WS_IC.disp_mod('w',cell2mat(cellfun(@(w)w(:),W_IC,'Un',0)));
Y_mod = WS_Yeq.disp_mod('w',cell2mat(cellfun(@(w)w(:),W_Y,'Un',0)));
X_mod = WS_Xeq.disp_mod('w',cell2mat(cellfun(@(w)w(:),W_X,'Un',0)));

disp(['IC terms:'])
IC_mod{:}
disp(['Yeq terms:'])
Y_mod{:}
disp(['Xeq terms:'])
X_mod{:}

%% view conf int
coeff_compare;

figure(1);clf
varget = 'IC';
view_conf_int;

figure(2);clf
varget = 'Y';
view_conf_int;

figure(3);clf
varget = 'X';
view_conf_int;

%% decide on terms to prune

IC_prune = [2 3 4 6];
Y_prune = [4 6 7];
X_prune = [1];

W_IC_sparse = zero_ents(W_IC,IC_prune);
W_Y_sparse = zero_ents(W_Y,Y_prune);
W_X_sparse = zero_ents(W_X,X_prune);

tags_IC_X = lib_X_IC.tags;
tags_Y_Yeq = lib_Y_Yeq.tags;
tags_X_Yeq = lib_X_Yeq.tags;
tags_X_Xeq = lib_X_Xeq.tags;
tags_Y_Xeq = lib_Y_Xeq.tags;

%% extract model terms from wsindy_eco run

WENDy_args = {'maxits',10,'ittol',10^-4,'diag_reg',10^-6,'verbose',1};
[rhs_IC,W_IC,rhs_Y,W_Y,rhs_X,W_X,...
    lib_Y_IC,lib_X_IC,lib_Y_Yeq,lib_X_Yeq,lib_Y_Xeq,lib_X_Xeq,...
    WS_IC,WS_Yeq,WS_Xeq,...
    CovW_IC,CovW_Y,CovW_X,...
    errs_Yend]= ...
    wendy_eco_fcn(toggle_zero_crossing,stop_tol,phifun_Y,tf_Y_params,WENDy_args,autowendy,tol_dd_learn,...
        Y_train,X_train,X_var,train_inds,train_time,t_epi,yearlength,nstates_X,nstates_Y,X_in,nX,nY,...
        tags_IC_X,tags_Y_Yeq,tags_X_Yeq,tags_X_Xeq,tags_Y_Xeq,W_IC_sparse,W_Y_sparse,W_X_sparse);

%% sim full system

sim_script;