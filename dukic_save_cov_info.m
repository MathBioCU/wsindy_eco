Cov_Res = WS_Yeq.cov;
w = WS_Yeq.weights;
w_sparse = w(w~=0);
[~,W_Y_c,~] = shorttime_map(W_Y_compare,lib_Y_Yeq,lib_X_Yeq,1./nX,1./nY);
w_true = cell2mat(cellfun(@(w)reshape(w',[],1),W_Y_c','uni',0));
w_true_sparse = w_true(w_true~=0);
res = res_0_Y(:,end);
res_true = WS_Yeq.G{1}*w_true-WS_Yeq.b{1};
G = WS_Yeq.G{1};
b = WS_Yeq.b{1};
CovW = CovW_Y;
save(['~/Desktop/WSINDy_eco_cov_info_TPR',num2str(tpr_Y),'.mat'],'Cov_Res','w','w_sparse','w_true','w_true_sparse','res','res_true','G','b','CovW')