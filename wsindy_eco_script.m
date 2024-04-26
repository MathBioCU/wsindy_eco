%% load hyperparams

wsindy_eco_inputs; tic;

%% run alg

[rhs_IC,W_IC,rhs_Y,W_Y,rhs_X,W_X,lib_Y_IC,lib_X_IC,...
    lib_Y_Yeq,lib_X_Yeq,lib_Y_Xeq,lib_X_Xeq,WS_IC,WS_Yeq,WS_Xeq,...
    CovW_IC,CovW_Y,CovW_X,errs_Yend]= ...
    wsindy_eco_fcn(toggle_zero_crossing,stop_tol,phifun_Y,tf_Y_params,WENDy_args,autowendy,tol,tol_min,tol_dd_learn,pmax_IC,polys_Y_Yeq,polys_X_Xeq,pmax_X_Yeq,pmax_Y_Xeq,neg_Y,neg_X,boolT,...
    Y_train,X_train,train_inds,train_time,t_epi,yearlength,nstates_X,nstates_Y,X_in,nX,nY,...
    custom_tags_X,custom_tags_Y,linregargs_fun_IC,linregargs_fun_Y,linregargs_fun_X);

%% coeff compare

W_IC_compare = inject_coeff_param(W_IC_true,zeros(1,nstates_Y),tags_IC_true,cell2mat(lib_Y_IC.tags'),cell2mat(lib_X_IC.tags'));
errs_2_IC = norm(reshape([W_IC{:}]-[W_IC_compare{:}],[],1))/norm(reshape([W_IC_compare{:}],[],1));
fprintf('IC coeff err 2: %1.3e \n',errs_2_IC);
errs_inf_IC = abs(reshape([W_IC{:}]-[W_IC_compare{:}],[],1))./abs(reshape([W_IC_compare{:}],[],1));
errs_inf_IC = max(errs_inf_IC(errs_inf_IC<inf),[],'all','omitnan');
if isempty(errs_inf_IC)
    errs_inf_IC = NaN;
end
fprintf('IC coeff errInf: %1.3e \n',errs_inf_IC)
tpr_IC = tpscore(reshape([W_IC{:}],[],1),reshape([W_IC_compare{:}],[],1));
fprintf('IC TPR: %0.3f \n',tpr_IC)

W_Y_compare = inject_coeff_param(W_Y_true,tags_Y_true,tags_X_true, cell2mat(lib_Y_Yeq.tags'),cell2mat(lib_X_Yeq.tags'));
errs_2_Y = norm(reshape([W_Y{:}]-[W_Y_compare{:}],[],1))/norm(reshape([W_Y_compare{:}],[],1));
fprintf('Ydot coeff err: %1.3e \n',errs_2_Y)
errs_inf_Y = abs(reshape([W_Y{:}]-[W_Y_compare{:}],[],1))./abs(reshape([W_Y_compare{:}],[],1));
errs_inf_Y = max(errs_inf_Y(errs_inf_Y<inf),[],'all','omitnan');
if isempty(errs_inf_Y)
    errs_inf_Y = NaN;
end
fprintf('Ydot coeff errInf: %1.3e \n',errs_inf_Y)
tpr_Y = tpscore(reshape([W_Y{:}],[],1),reshape([W_Y_compare{:}],[],1));
fprintf('Ydot TPR: %0.3f \n',tpr_Y)

W_X_compare = inject_coeff_param(W_X_true,tags_Ext_X_true,tags_Ext_Y_true,cell2mat(lib_X_Xeq.tags'),cell2mat(lib_Y_Xeq.tags'));
errs_2_X = norm(reshape([W_X{:}]-[W_X_compare{:}],[],1))/norm(reshape([W_X_compare{:}],[],1));
fprintf('Xn coeff err2: %1.3e \n',errs_2_X)
errs_inf_X = abs(reshape([W_X{:}]-[W_X_compare{:}],[],1))./abs(reshape([W_X_compare{:}],[],1));
errs_inf_X = max(errs_inf_X(errs_inf_X<inf),[],'all','omitnan');
if isempty(errs_inf_X)
    errs_inf_X = NaN;
end
fprintf('Xn coeff errInf: %1.3e \n',errs_inf_X)
tpr_X = tpscore(reshape([W_X{:}],[],1),reshape([W_X_compare{:}],[],1));
fprintf('Xn TPR: %0.3f \n',tpr_X)
fprintf('\n runtime: %2.3f \n',toc)

%% sim full system
addpath(genpath('wsindy_obj_base'))
if toggle_sim>0
    x0s = [X(1,:);mean(X_train.*nX).*(1 + sqrt(3)*0.2*(rand(num_sim,2)-0.5)*2)];
    for j=1:size(x0s,1)
        x0 = x0s(j,:);
        num_gen = floor(size(X,1));
        sig_tmax = 0;
        num_t_epi = 112;
        if isequal(x0,X(1,:))
            X_test = X; Ycell_test = Ycell; Y_test = Y; t_epi_test = t_epi; tn_test = tn; t_test = t;
        else
            [X_test,Ycell_test,Y_test,t_test,tn_test,t_epi_test] = sim_hybrid_fcn(rhs_IC_true,rhs_Y_true,rhs_X_true,x0,nstates_Y,...
            num_gen,num_t_epi,yearlength,sig_tmax,tol_dd_sim,toggle_zero_crossing,inf);
        end
        [X_pred,Ycell_pred,Y_pred,t_pred,tn_pred,t_epi_pred] = sim_hybrid_fcn(rhs_IC,rhs_Y,rhs_X,x0,nstates_Y,...
            num_gen,num_t_epi,yearlength,sig_tmax,tol_dd_sim,toggle_zero_crossing,stop_tol*max(max(cell2mat(Y_train)./nY)));
    
        n = size(X_pred,1);
        if toggle_zero_crossing==1
            if any([X_pred(n,:)]<0)
                fprintf(['negative values detected in X(%u) \n'],find([X_pred(n,:)]<0))
            end
            if any(Ycell_pred{n-1}(:)<0)
                fprintf(['negative values detected in Y, %u \n'],find(Ycell_pred{n-1}(:)<0))
            end
        end
        
        fprintf('X err: %1.3e \n',vecnorm(X_pred(1:n,:)-X_test(1:n,:))./vecnorm(X_test(1:n,:)))
        cumerr = arrayfun(@(i)norm(vecnorm(X_pred(1:i,:)-X_test(1:i,:),2,2))/norm(vecnorm(X_test(1:i,:),2,2)),(1:n)');
        n_err_tol = find(cumerr>err_tol,1);
        if isempty(n_err_tol)
            n_err_tol = n-1;
        else
            n_err_tol = n_err_tol-1;
        end
        fprintf('n_{%0.1f}=%0.2f \n',err_tol,n_err_tol)
        
        if toggle_vis
            figure(j);clf;set(gcf,'position',[1755         156        1837         716])
            wsindy_eco_vis;
        end
    end
end
%% functions

function [value, isterminal, direction] = myEvent(T, Y, thresh, toggle_zero_crossing)
    if toggle_zero_crossing
        value      = or(norm(Y) >= thresh, any(Y==0));
    else
        value      = norm(Y) >= thresh;
    end
    isterminal = 1;   % Stop the integration
    direction  = 0;
end