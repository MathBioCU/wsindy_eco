x_cl = 'ko'; % colors for true data
y_cl = 'r'; % colors for true data
xL_cl = 'bo'; % colors for learned data
yL_cl = 'c--';% colors for learned data
xO_cl = 'ko';% colors for observed data
yO_cl = 'r'; % colors for observed data
ylims = arrayfun(@(j)[0.1*min([X_pred(:,j);X_test(:,j)]) 10*max([X_pred(:,j);X_test(:,j)])],1:nstates_X,'uni',0);
yticks = arrayfun(@(j)10.^(floor(log10(ylims{j}(1))):2:ceil(log10(ylims{j}(2)))),1:nstates_X,'uni',0);
for j=1:nstates_X
    subplot(nstates_X,2,(j-1)*nstates_X+1)
    plot(tn_test(1:n),X_test(1:n,j),x_cl,tn_pred(1:n),X_pred(1:n,j),xL_cl,...
        t_test, Y_test(:,j),y_cl,t_pred, Y_pred(:,j),yL_cl,'linewidth',3,'markersize',7)
    set(gca,'ticklabelinterpreter','latex','fontsize',20,...
        'Xtick',yearlength*floor(linspace(0,max(n-1,1),min(n,5))),...
        'XtickLabels',floor(linspace(0,max(n-1,1),min(n,5))),...
        'Yscale','log','Ylim',ylims{j},'Ytick',yticks{j},'Xlim',[0 yearlength*max(n-1,1)])
    grid on;
    if j==1
        legend({'${\bf X_test^\star}$','$\widehat{\bf X_test}$','${\bf Y_test^\star}$','$\widehat{\bf Y_test}$'},'location','nw','interpreter','latex','fontsize',16)
    end
    ylabel(['Compartment ',num2str(j)],'Rotation',90,'interpreter','latex')
    if j==nstates_X
        xlabel('$n$','interpreter','latex')
    end
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
                hold off
            end
        end    
        if j~=1
            hold on
            yl = max(get(gca,'ylim'),10^-5);
            h3 = plot(yearlength*(n_err_tol*[1 1]-1),[yl(1) yl(2)],'k--','linewidth',3);
            hold off
            legend([h1;h2;h3],{'${\bf X_test}$','${\bf Y_test}$','$n_{0.2}(\hat{\bf w})$'},'location','ne','interpreter','latex','fontsize',16)
        end
        if j==nstates_X
            xlabel('$n$','interpreter','latex')        
        end
        set(gca,'ticklabelinterpreter','latex','fontsize',20,...
            'Xtick',yearlength*floor(linspace(0,max(n-1,1),min(n,5))),...
            'XtickLabels',floor(linspace(0,max(n-1,1),min(n,5))),...
            'Yscale','log','Ylim',ylims{j},'Ytick',yticks{j},'Xlim',[0 yearlength*max(n-1,1)])
        grid on
        % if j==1
        %     [~,a]=sort(cellfun(@(Y_test)max(Y_test(:,j)),Y_train),'descend');
        %     axes('Position',[.76 0.65 .13 .24])
        %     box on    
        %     hold on
        %     ks = a(5);
        %     for k=ks
        %         i = train_inds(X_in(k));
        %         plot((i-1)*yearlength+t_epi_test{i},...
        %             Ycell_test{i}(:,j),yO_cl,'markersize',7,'linewidth',4)
        %         plot((i-1)*yearlength+train_time{i==train_inds(X_in)},...
        %             Y_train{i==train_inds(X_in)}(:,j)*nY(j),'m-','markersize',7,'linewidth',4)
        %         legend({'${\bf Y_test^\star}$','${\bf Y_test}$'},'interpreter','latex','fontsize',16)
        %     end
        %     xlim([(train_inds(X_in(min(ks)))-1)*yearlength train_inds(X_in(max(ks)))*yearlength])
        %     hold off
        %     set(gca,'ticklabelinterpreter','latex','fontsize',16,...
        %         'Xtick',yearlength*linspace(0,max(n-1,1),min(n,5)*50),...
        %         'XtickLabels',round(linspace(0,max(n-1,1),min(n,5)*50),2),...
        %         'Yscale','log')
        %     grid on;
        % end
    end
else
    for j=1:nstates_X
        subplot(nstates_X,2,(j-1)*nstates_X+2)
        plot(tn_test,X_test(:,j),x_cl,t_test, Y_test(:,j),y_cl,'linewidth',3,'markersize',7)
        xlim([0 yearlength*(size(X_test,1)-1)])
        set(gca,'ticklabelinterpreter','latex','fontsize',16,...
                'Yscale','log','Ylim',ylims{j},'Ytick',yticks{j})
        grid on;
        if j==nstates_X
            xlabel('$n$','interpreter','latex')        
        end
    end
end
saveas(gcf,['~/Desktop/WENDy_hybrid_sol_snr',num2str(snr_Y),'.png'])