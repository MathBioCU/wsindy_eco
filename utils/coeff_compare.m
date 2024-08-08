%% coeff compare

W_IC_compare = inject_coeff_param(W_IC_true(:),zeros(1,nstates_Y),tags_IC_true,cell2mat(lib_Y_IC.tags'),cell2mat(lib_X_IC.tags'));
errs_2_IC = cellfun(@(W,V) norm(W-V)/norm(V),W_IC,W_IC_compare);
errs_inf_IC = cellfun(@(W,V) max(abs(W(V~=0)-V(V~=0))./abs(V(V~=0))),W_IC,W_IC_compare,'un',0);
errs_inf_IC = cell2mat(errs_inf_IC(cellfun(@(e)~isempty(e),errs_inf_IC)));
tpr_IC = cellfun(@(W,V)tpscore(W,V),W_IC,W_IC_compare);
fprintf('IC coeff err 2: %1.3e \n',errs_2_IC);
fprintf('IC coeff errInf: %1.3e \n',errs_inf_IC)
fprintf('IC TPR: %0.3f \n',tpr_IC)

if length(lib_Y_Yeq)==1
    lib_Y_Yeq = repmat(lib_Y_Yeq,1,nstates_Y);
end
W_Y_compare = arrayfun(@(i)...
    cell2mat(inject_coeff_param(W_Y_true(i),tags_Y_true,tags_X_true, cell2mat(lib_Y_Yeq(i).tags'),cell2mat(lib_X_Yeq.tags'))),...
    (1:length(W_Y_true))','un',0);
errs_2_Y = cellfun(@(W,V) norm(W-V)/norm(V),W_Y,W_Y_compare);
errs_inf_Y = cellfun(@(W,V) max(abs(W(V~=0)-V(V~=0))./abs(V(V~=0))),W_Y,W_Y_compare,'un',0);
errs_inf_Y = cell2mat(errs_inf_Y(cellfun(@(e)~isempty(e),errs_inf_Y)));
tpr_Y = cellfun(@(W,V)tpscore(W,V),W_Y,W_Y_compare);
fprintf('Ydot coeff err: %1.3e \n',errs_2_Y)
fprintf('Ydot coeff errInf: %1.3e \n',errs_inf_Y)
fprintf('Ydot TPR: %0.3f \n',tpr_Y)

W_X_compare = inject_coeff_param(W_X_true(:),tags_Ext_X_true,tags_Ext_Y_true,cell2mat(lib_X_Xeq.tags'),cell2mat(lib_Y_Xeq.tags'));
errs_2_X = cellfun(@(W,V) norm(W-V)/norm(V),W_X,W_X_compare);
errs_inf_X = cellfun(@(W,V) max(abs(W(V~=0)-V(V~=0))./abs(V(V~=0))),W_X,W_X_compare,'un',0);
errs_inf_X = cell2mat(errs_inf_X(cellfun(@(e)~isempty(e),errs_inf_X)));
tpr_X = cellfun(@(W,V)tpscore(W,V),W_X,W_X_compare);
fprintf('Xn coeff err2: %1.3e \n',errs_2_X)
fprintf('Xn coeff errInf: %1.3e \n',errs_inf_X)
fprintf('Xn TPR: %0.3f \n',tpr_X)
