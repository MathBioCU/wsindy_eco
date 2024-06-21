
function W_s = wendy_hybrid_sample(W,C)
    [W_p,ss] = wendy_param(W);
    W_samp = mvnrnd(W_p,C);
    W_s = W;
    ind = 0;
    for j=1:length(W)
        W_s{j}(ss{j}) = W_samp(ind+1:ind+length(ss{j}));
        ind = ind + length(ss{j});
    end
end
