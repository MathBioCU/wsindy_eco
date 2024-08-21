addpath(genpath('wsindy_obj_base'))
if ~exist('yscl','var')
    yscl = 'log';
end
if toggle_sim>0
    %%% get initial conditions to test on
    x0s = [X(train_inds(1),:);X(train_inds(end)+1,:);mean(X_train.*nX).*(1 + sqrt(3)*oos_std*(rand(num_sim,nstates_X)-0.5)*2)];

    %%% save test sims
    t_test_cell = cell(size(x0s,1),1);
    X_test_cell = cell(size(x0s,1),1);
    Y_test_cell = cell(size(x0s,1),1);
    t_pred_cell = cell(size(x0s,1),1);
    X_pred_cell = cell(size(x0s,1),1);
    Y_pred_cell = cell(size(x0s,1),1);

    %%% simulate learned model and ground truth if available
    for jj=1:size(x0s,1)
        x0 = x0s(jj,:); % current IC
        num_t_epi = size(Ycell{1},1); % number of timepoints

        %%% get ground truth data
        if isequal(x0,X(1,:)) % training set
            X_test_cell{jj} = X; 
            Y_test_cell{jj} = Ycell; 
            Y_test = Y; 
            t_epi_test = t_epi; 
            tn_test = tn; 
            t_test_cell{jj} = t;
        elseif isequal(x0,X(train_inds(end)+1,:)) % testing set 
            X_test_cell{jj}=X(train_inds(end)+1:end,:);
            Y_test_cell{jj}=Ycell(train_inds(end)+1:end);
            t_epi_test = t_epi(train_inds(end)+1:end);
            tn_test = (0:num_gen-train_inds(end)-1)*yearlength;
            t_test = cell2mat(arrayfun(@(i)(i-1)*yearlength+t_epi_test{i},(1:length(t_epi_test))','uni',0)); % full continuous time
            Y_test = cell2mat( Y_test_cell{jj});
            t_test_cell{jj} = t_test;
        else % validation test on random IC
            if exist('rhs_IC_true','var') % if black box simulator available, simulate new ground truth data
                [X_test_cell{jj},Y_test_cell{jj},Y_test,t_test_cell{jj},tn_test,t_epi_test] = sim_hybrid_fcn(...
                        rhs_IC_true,rhs_Y_true,rhs_X_true,x0,nstates_Y,...
                        num_gen,num_t_epi,yearlength,0,tol_dd_sim,toggle_zero_crossing,inf);

            else 
                X_test_cell{jj} = [];
                Y_test_cell{jj} = [];
                t_test_cell{jj} = [];
            end
        end

        %%% get learned model output
        [X_pred_cell{jj},Y_pred_cell{jj},Y_pred,t_pred_cell{jj},tn_pred,t_epi_pred] = sim_hybrid_fcn(...
                rhs_IC,rhs_Y,rhs_X,x0,nstates_Y,...
                num_gen,num_t_epi,yearlength,0,tol_dd_sim,toggle_zero_crossing,stop_tol*max(max(cell2mat(Y_train)./nY)));
    
        %%% assess performance
        n = min(size(X_pred_cell{jj},1),size(X_test_cell{jj},1));
        n_err_tol = get_n_err_tol(X_pred_cell{jj},X_test_cell{jj},err_tol);
        if n>0
            fprintf('X err: %1.3e \n',vecnorm(X_pred_cell{jj}(1:n,:)-X_test_cell{jj}(1:n,:))./vecnorm(X_test_cell{jj}(1:n,:)))
            fprintf('n_{%0.1f}=%0.2f \n',err_tol,n_err_tol)
        end
        
        %%% detect zero crossings, get 
        if toggle_zero_crossing==1
            if any(X_pred_cell{jj}<0)
                gennum = find(sum(X_pred_cell{jj}<0,2)==1,1);
                fprintf(['negative values detected in X(%u) \n'],X_pred_cell{jj}(gennum,:))
            end
            if any(Y_pred_cell{jj}{end}(:)<0)
                gennum = find(cellfun(@(y) any(sum(y<0,2)==1),Y_pred_cell{jj}),1);
                fprintf(['negative values detected in Y(%u) \n'],find(sum(Y_pred_cell{jj}{gennum}<0)))
            end
        end

        %%% get yticks and ylims
        if toggle_zero_crossing==1
            if exist('rhs_IC_true','var')
                ylims = arrayfun(@(j)[0.1*min([X_pred_cell{jj}(:,j);X_test_cell{jj}(:,j)]) 10*max([X_pred_cell{jj}(:,j);X_test_cell{jj}(:,j)])],...
                    1:nstates_X,'uni',0);
            else
                ylims = arrayfun(@(j)[0.1*min([X_pred_cell{jj}(1:end-1,j)]) 10*max([X_pred_cell{jj}(1:end-1,j)])],1:nstates_X,'uni',0);
            end
            yticks = arrayfun(@(j)10.^(floor(log10(max(ylims{j}(1),eps))):2:ceil(log10(ylims{j}(2)))),1:nstates_X,'uni',0);
        else
            ylims = arrayfun(@(j)[-inf inf],1:nstates_X,'uni',0);
            yticks = arrayfun(@(j)[],1:nstates_X,'uni',0);
        end
        
        if toggle_vis==1
            figure(jj);clf;
            x_cl = 'ko'; % colors for true data
            y_cl = 1.3*[0.4660 0.6740 0.1880]; % colors for true data
            xL_cl = 'ko'; % colors for learned data
            yL_cl = 'b-.';% colors for learned data
            xO_cl = 'ko';% colors for observed data
            yO_cl = 'r'; % colors for observed data
            yobs_cl = 'r-';
            
            for j=1:nstates_X
                subplot(1,nstates_X,j)
                    %%% plot ground truth data
                    if n>0
                        h0 = plot(tn_test,X_test_cell{jj}(:,j),x_cl,...
                                  t_test_cell{jj}, Y_test(:,j),'-','linewidth',3,'markersize',7);
                        set(h0(2),'color',y_cl)
                        hold on
                        plot(yearlength*[n_err_tol]*[1 1],ylims{j},'k--','linewidth',3)
                    end
    
                    h2 = plot(tn_pred,X_pred_cell{jj}(:,j),xL_cl,...
                              t_pred_cell{jj}, Y_pred(:,j),yL_cl,'linewidth',3,'markersize',7);

                    set(gca,'ticklabelinterpreter','latex','fontsize',12,...
                        'Xtick',yearlength*(0:ceil(num_gen/6):num_gen),...
                        'XtickLabels',0:ceil(num_gen/6):num_gen,...
                        'Yscale',yscl,...
                        'Ylim',ylims{j},...
                        'Ytick',yticks{j});%,...
                        ...'Xlim',[0 yearlength*max(n-2,1)]
                    grid on;
                
                    if exist('h0','var')
                        legend([h0(2);h2(2)],{'true model output','learned model output'},'location','sw','interpreter','latex','fontsize',14)
                    else
                        legend(h2(2),{'learned model output'},'location','sw','interpreter','latex','fontsize',14)
                    end

                    ylabel(['X(',num2str(j),')'])
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