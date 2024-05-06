addpath(genpath('wsindy_obj_base'))
if toggle_sim>0
    x0s = [X(train_inds(1),:);mean(X_train.*nX).*(1 + sqrt(3)*oos_std*(rand(num_sim,2)-0.5)*2)];
    for jj=1:size(x0s,1)
        x0 = x0s(jj,:);
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
                fprintf(['negative values detected in Y(%u) \n'],find(Ycell_pred{n-1}(end,:)<0))
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
            figure(jj);clf;
            % set(gcf,'position',[1755         156        1837         716])
            wsindy_eco_vis;
        end
    end
end