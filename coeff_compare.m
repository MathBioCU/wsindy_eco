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
