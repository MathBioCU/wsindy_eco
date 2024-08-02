W_IC_compare = inject_coeff_param(W_IC_true(:),zeros(1,nstates_Y),tags_IC_true,cell2mat(lib_Y_IC.tags'),cell2mat(lib_X_IC.tags'));
W_Y_compare = arrayfun(@(i)...
    cell2mat(inject_coeff_param(W_Y_true(i),tags_Y_true,tags_X_true, cell2mat(lib_Y_Yeq(i).tags'),cell2mat(lib_X_Yeq.tags'))),...
    (1:length(W_Y_true))','un',0);
W_X_compare = inject_coeff_param(W_X_true(:),tags_Ext_X_true,tags_Ext_Y_true,cell2mat(lib_X_Xeq.tags'),cell2mat(lib_Y_Xeq.tags'));



W_Y_compvec = cell2mat(cellfun(@(w)w(:),W_Y_compare,'uni',0));
S = {W_Y_compvec~=0};
WS_Yeq.weights = [];
WS = WS_opt().wendy(WS_Yeq,'linregargs',{'S',S});
W = arrayfun(@(i,L)reshape(WS.reshape_w{i},length(L.terms),[]),(1:WS.numeq)',lib_Y_Yeq','uni',0);
[rhs,W,~,M] = shorttime_map(W,lib_Y_Yeq,lib_X_Yeq,nY,nX);