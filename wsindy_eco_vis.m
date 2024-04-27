x_cl = 'ko'; % colors for true data
y_cl = 'r'; % colors for true data
xL_cl = 'bo'; % colors for learned data
yL_cl = 'c--';% colors for learned data
xO_cl = 'ko';% colors for observed data
yO_cl = 'r'; % colors for observed data
if toggle_zero_crossing
    ylims = arrayfun(@(j)[0.1*min([X_pred(:,j);X_test(:,j)]) 10*max([X_pred(:,j);X_test(:,j)])],1:nstates_X,'uni',0);
    yticks = arrayfun(@(j)10.^(floor(log10(ylims{j}(1))):2:ceil(log10(ylims{j}(2)))),1:nstates_X,'uni',0);
    YS = 'log';
else
    ylims = arrayfun(@(j)[-inf inf],1:nstates_X,'uni',0);
    yticks = arrayfun(@(j)[],1:nstates_X,'uni',0);
    YS = 'linear';
end
for j=1:nstates_X
    subplot(nstates_X,2,(j-1)*nstates_X+1)
    hold on
    plot(tn_test(1:n),X_test(1:n,j),x_cl,tn_pred(1:n),X_pred(1:n,j),xL_cl,...
        t_test, Y_test(:,j),y_cl,t_pred, Y_pred(:,j),yL_cl,'linewidth',3,'markersize',7)
    if toggle_zero_crossing
        yl = max(ylims{j},min(X_pred(X_pred>0)));
    else
        yl = get(gca,'Ylim');
    end
    h3 = plot(yearlength*(n_err_tol*[1 1]-1),yl,'k--','linewidth',3);
    hold off
    set(gca,'ticklabelinterpreter','latex','fontsize',20,...
        'Xtick',yearlength*floor(linspace(0,max(n-1,1),min(n,5))),...
        'XtickLabels',floor(linspace(0,max(n-1,1),min(n,5))),...
        'Yscale',YS,'Ylim',ylims{j},'Ytick',yticks{j},'Xlim',[0 yearlength*max(n-1,1)])
    grid on;
    legend({'$X_{test}$','$X_{pred}$','$Y_{test}$','$Y_{pred}$','$n_{0.2}(\hat{\bf w})$'},'location','bestoutside','interpreter','latex','fontsize',18)
    ylabel(['Compartment ',num2str(j)],'Rotation',90,'interpreter','latex')
    xlabel('$n$','interpreter','latex')
end
if isequal(x0,X(1,:))
    for j=1:nstates_X
        for i=1:size(X_test,1)
            if ismember(i,train_inds)
                subplot(nstates_X,2,(j-1)*nstates_X+2)
                hold on
                h1=plot((i-1)*yearlength,X_train(i==train_inds,j)*nX(j),xO_cl,'markersize',7,'linewidth',5);
                if ismember(i,train_inds(X_in))
                    h2=plot((i-1)*yearlength+train_time{i==train_inds(X_in)},...
                        Y_train{i==train_inds(X_in)}(:,j)*nY(j),'m-','markersize',7,'linewidth',4);
                end
            end
        end    
        hold on
        hold off
        legend([h1;h2],{'$X_{train}$','$Y_{train}$','$n_{0.2}(\hat{\bf w})$'},'location','ne','interpreter','latex','fontsize',16)
        xlabel('$n$','interpreter','latex')        
        set(gca,'ticklabelinterpreter','latex','fontsize',20,...
            'Xtick',yearlength*floor(linspace(0,max(n-1,1),min(n,5))),...
            'XtickLabels',floor(linspace(0,max(n-1,1),min(n,5))),...
            'Yscale',YS,'Ylim',ylims{j},'Ytick',yticks{j},'Xlim',[0 yearlength*max(n-1,1)])
        grid on
    end
else
    for j=1:nstates_X
        subplot(nstates_X,2,(j-1)*nstates_X+2)
        plot(tn_test,X_test(:,j),x_cl,t_test, Y_test(:,j),y_cl,'linewidth',3,'markersize',7)
        xlim([0 yearlength*(size(X_test,1)-1)])
        set(gca,'ticklabelinterpreter','latex','fontsize',16,...
                'Yscale',YS,'Ylim',ylims{j},'Ytick',yticks{j})
        grid on;
        xlabel('$n$','interpreter','latex')        
    end
end
% saveas(gcf,['~/Desktop/WENDy_hybrid_sol_snr',num2str(snr_Y),'.png'])