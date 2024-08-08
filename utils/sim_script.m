addpath(genpath('wsindy_obj_base'))
if toggle_sim>0
    x0s = [X(train_inds(1),:);mean(X_train.*nX).*(1 + sqrt(3)*oos_std*(rand(num_sim,nstates_X)-0.5)*2)];
    t_test_cell = cell(size(x0s,1),1);
    X_test_cell = cell(size(x0s,1),1);
    Y_test_cell = cell(size(x0s,1),1);
    t_pred_cell = cell(size(x0s,1),1);
    X_pred_cell = cell(size(x0s,1),1);
    Y_pred_cell = cell(size(x0s,1),1);
    for jj=1:size(x0s,1)
        x0 = x0s(jj,:);
        num_gen = floor(size(X,1));
        sig_tmax = 0;
        num_t_epi = size(Ycell{1},1);
        if isequal(x0,X(1,:))
            X_test_cell{jj} = X; Y_test_cell{jj} = Ycell; Y_test = Y; t_epi_test = t_epi; tn_test = tn; t_test_cell{jj} = t;
        else
            if exist('rhs_IC_true','var')
                [X_test_cell{jj},Y_test_cell{jj},Y_test,t_test_cell{jj},tn_test,t_epi_test] = sim_hybrid_fcn(rhs_IC_true,rhs_Y_true,rhs_X_true,x0,nstates_Y,...
                num_gen,num_t_epi,yearlength,sig_tmax,tol_dd_sim,toggle_zero_crossing,inf);
            else
                X_test_cell{jj}=X*0;
                Y_test_cell{jj}=cellfun(@(y)y*0,Ycell,'uni',0);
                t_epi_test = t_epi; tn_test = tn; t_test = t;
                Y_test = Y*0;
            end
        end
        [X_pred_cell{jj},Y_pred_cell{jj},Y_pred,t_pred_cell{jj},tn_pred,t_epi_pred] = sim_hybrid_fcn(rhs_IC,rhs_Y,rhs_X,x0,nstates_Y,...
            num_gen,num_t_epi,yearlength,sig_tmax,tol_dd_sim,toggle_zero_crossing,stop_tol*max(max(cell2mat(Y_train)./nY)));
    
        n = min(size(X_pred_cell{jj},1),size(X_test_cell{jj},1));
        if toggle_zero_crossing==1
            if any([X_pred_cell{jj}(n,:)]<0)
                fprintf(['negative values detected in X(%u) \n'],find([X_pred_cell{jj}(n,:)]<0))
            end
            if any(Y_pred_cell{jj}{n-1}(:)<0)
                fprintf(['negative values detected in Y(%u) \n'],find(Y_pred_cell{jj}{n-1}(end,:)<0))
            end
        end
        
        fprintf('X err: %1.3e \n',vecnorm(X_pred_cell{jj}(1:n,:)-X_test_cell{jj}(1:n,:))./vecnorm(X_test_cell{jj}(1:n,:)))
        cumerr = arrayfun(@(i)norm(vecnorm(X_pred_cell{jj}(1:i,:)-X_test_cell{jj}(1:i,:),2,2))/norm(vecnorm(X_test_cell{jj}(1:i,:),2,2)),(1:n)');
        n_err_tol = get_n_err_tol(X_pred_cell{jj},X_test_cell{jj},err_tol);
        fprintf('n_{%0.1f}=%0.2f \n',err_tol,n_err_tol)

        if toggle_zero_crossing
            if exist('rhs_IC_true','var')
                ylims = arrayfun(@(j)[0.1*min([X_pred_cell{jj}(:,j);X_test_cell{jj}(:,j)]) 10*max([X_pred_cell{jj}(:,j);X_test_cell{jj}(:,j)])],1:nstates_X,'uni',0);
            else
                ylims = arrayfun(@(j)[0.1*min([X_pred_cell{jj}(1:end-1,j)]) 10*max([X_pred_cell{jj}(1:end-1,j)])],1:nstates_X,'uni',0);
            end
            yticks = arrayfun(@(j)10.^(floor(log10(max(ylims{j}(1),eps))):2:ceil(log10(ylims{j}(2)))),1:nstates_X,'uni',0);
            YS = 'log';
        else
            ylims = arrayfun(@(j)[-inf inf],1:nstates_X,'uni',0);
            yticks = arrayfun(@(j)[],1:nstates_X,'uni',0);
            YS = 'linear';
        end
        
        if toggle_vis
            figure(jj);clf;
            x_cl = 'ko'; % colors for true data
            y_cl = 1.3*[0.4660 0.6740 0.1880]; % colors for true data
            xL_cl = 'ko'; % colors for learned data
            yL_cl = 'b--';% colors for learned data
            xO_cl = 'ko';% colors for observed data
            yO_cl = 'r'; % colors for observed data
            yobs_cl = 'r-';
            
            for j=1:nstates_X
                subplot(1,nstates_X,j)
                h0=plot(tn_test(1:n),X_test_cell{jj}(1:n,j),x_cl,t_test_cell{jj}, Y_test(:,j),'-','linewidth',3,'markersize',7);
                hold on
                if isequal(x0,X(1,:))
                    plot(yearlength*[n_err_tol]*[1 1],ylims{j},'k--','linewidth',3)
                end
                set(h0(2),'color',y_cl)
                h2=plot(tn_pred(1:n),X_pred_cell{jj}(1:n,j),xL_cl,t_pred_cell{jj}, Y_pred(:,j),yL_cl,'linewidth',3,'markersize',7);
                if isequal(x0,X(1,:))
                    plot(yearlength*[n_err_tol]*[1 1],ylims{j},'k--','linewidth',3)
                end                
                set(gca,'ticklabelinterpreter','latex','fontsize',12,...
                    'Xtick',yearlength*floor(linspace(0,max(n-1,1),min(n,5))),...
                    'XtickLabels',floor(linspace(0,max(n-1,1),min(n,5))),...
                    'Yscale',yscl,...
                    ...'Ylim',ylims{j},'Ytick',yticks{j},...
                    'Xlim',[0 yearlength*max(n-2,1)])
                grid on;
            
                legend([h0(2);h2(2)],{'true model output','learned model output'},'location','sw','interpreter','latex','fontsize',14)
                if j==1
                    ylabel(['host $(S,N)$'],'interpreter','latex')
                else
                    ylabel(['pathogen $(P,Z)$'],'interpreter','latex')
                end
                xlabel('generation number ($n$)','interpreter','latex')   
            end
                
        end
    end
end