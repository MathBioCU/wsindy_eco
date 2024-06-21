
function [W_p,ss] = wendy_param(W)
    W_p = cell2mat(cellfun(@(w)w(:),W(:),'Un',0));
    W_p = W_p(W_p~=0);
    ss = cellfun(@(w)find(w),W,'Un',0);
end