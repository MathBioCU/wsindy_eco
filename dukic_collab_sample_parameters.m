%% get samples
n=96;
X_pred_cloud = cell(n,1);
x0 = X(1,:);    
parfor j=1:n
    rng('shuffle')
    W_IC_s = wendy_hybrid_sample(W_IC,CovW_IC);
    W_Y_s = wendy_hybrid_sample(W_Y,CovW_Y);
    W_X_s = wendy_hybrid_sample(W_X,CovW_X);
    
    rhs_IC_s = shorttime_map(W_IC_s,lib_Y_IC,lib_X_IC,ones(1,nstates_Y),ones(1,nstates_X));
    rhs_IC_s = @(X) rhs_IC_s(zeros(nstates_Y,1),X(:));
    rhs_Y_s = shorttime_map(W_Y_s,lib_Y_Yeq,lib_X_Yeq,ones(1,nstates_Y),ones(1,nstates_X));
    rhs_X_s = longtime_map(W_X_s,lib_X_Xeq,lib_Y_Xeq,ones(1,nstates_X),ones(1,nstates_Y));
    
    [X_pred_cloud{j},Ycell_pred,Y_pred,t_pred,tn_pred,t_epi_pred] = sim_hybrid_fcn(...
        rhs_IC_s,rhs_Y_s,rhs_X_s,x0,nstates_Y,...
        num_gen,num_t_epi,yearlength,sig_tmax,...
        tol_dd_sim,toggle_zero_crossing,...
        stop_tol*max(max(cell2mat(Y_train)./nY)));    
end

%% view 
figure(1);clf
hold on
cellfun(@(x)semilogy(x(:,1),'b'),X_pred_cloud)
semilogy(X(:,1),'linewidth',3)
set(gca,'Yscale','log')
hold off

figure(2);clf
hold on
cellfun(@(x)semilogy(x(:,2),'b'),X_pred_cloud)
semilogy(X(:,2),'linewidth',3)
set(gca,'Yscale','log')
hold off

function W_s = wendy_hybrid_sample(W,C)
    W_p = cell2mat(cellfun(@(w)w(:),W,'Un',0));
    W_p = W_p(W_p~=0);
    W_samp = mvnrnd(W_p,C);
    ss = cellfun(@(w)find(w),W,'Un',0);
    W_s = W;
    ind = 0;
    for j=1:length(W)
        W_s{j}(ss{j}) = W_samp(ind+1:ind+length(ss{j}));
        ind = ind + length(ss{j});
    end
end
