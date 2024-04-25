% close all;
x_cl = 'ko'; % colors for true data
y_cl = [0.4660 0.6740 0.1880]; % colors for true data
xL_cl = 'ko'; % colors for learned data
yL_cl = 'b-';% colors for learned data
xO_cl = 'ko';% colors for observed data
yO_cl = 'r'; % colors for observed data
yobs_cl = 'r-';
yscl = 'log';
% ylims = arrayfun(@(j)[0.01*min([X_pred(:,j);X_test(:,j)]) 10*max([X_pred(:,j);X_test(:,j)])],1:nstates_X,'uni',0);
if isequal(yscl,'log')
    ylims = arrayfun(@(j)[0.01*min(X_test(1:n,j)) 2*max(X_test(1:n,j))],1:nstates_X,'uni',0);
else
    ylims = arrayfun(@(j)[0.1*min(X_test(1:n,j)) max(X_test(1:n,j))],1:nstates_X,'uni',0);
end
yticks = arrayfun(@(j)10.^(floor(log10(ylims{j}(1))):ceil(log10(ylims{j}(2)))),1:nstates_X,'uni',0);
for j=1:nstates_X
    % subplot(nstates_X,2,(j-1)*nstates_X+1)
    figure((j-1)*nstates_X+1)
    h0=plot(tn_test(1:n),X_test(1:n,j),x_cl, ...
        t_test, Y_test(:,j),yearlength*[n_err_tol]*[1 1],ylims{j},'k--','linewidth',3,'markersize',7);
     set(gca,'ticklabelinterpreter','latex','fontsize',20,...
        'Xtick',yearlength*floor(linspace(0,max(n-1,1),min(n,5))),...
        'XtickLabels',floor(linspace(0,max(n-1,1),min(n,5))),...
        'Yscale',yscl,'Ylim',ylims{j},'Ytick',yticks{j},'Xlim',[0 yearlength*max(n-1,1)])
     set(h0(2),'color',y_cl)
     grid on
    % if j==1
    % end
    xlabel('generation number ($n$)','interpreter','latex')
    for i=1:size(X_test,1)
        if ismember(i,train_inds)
            % subplot(nstates_X,2,(j-1)*nstates_X+2)
            hold on
            h1=plot((i-1)*yearlength,X_train(i==train_inds,j)*nX(j),xO_cl,'markersize',7,'linewidth',5);
            if ismember(i,train_inds(X_in))
                h2=plot((i-1)*yearlength+train_time{i==train_inds(X_in)},...
                    Y_train{i==train_inds(X_in)}(:,j)*nY(j),yobs_cl,'markersize',7,'linewidth',4);
            end
            hold off
        end
    end    
    if j==1
        ylabel(['host $(S,N)$'],'interpreter','latex')
    else
        ylabel(['pathogen $(P,Z)$'],'interpreter','latex')
    end

    if isequal(yscl,'log')
        legend([h0(2);h2],{'true data','observed data'},'location','sw','interpreter','latex','fontsize',16)
    else
        legend([h0(2);h2],{'true data','observed data'},'location','nw','interpreter','latex','fontsize',16)
    end
    saveas(gcf,['~/Desktop/hybrid_data',num2str(j),'_snr',num2str(snr_Y),'_numI',num2str(num_train_inds),'_subt,',num2str(subsamp_t),'_rng',num2str(seed1),'.png'])

    figure((j-1)*nstates_X+2)
    h2=plot(tn_pred(1:n),X_pred(1:n,j),xL_cl,...
        t_pred, Y_pred(:,j),yL_cl,yearlength*[n_err_tol]*[1 1],ylims{j},'k--','linewidth',3,'markersize',7);
    set(gca,'ticklabelinterpreter','latex','fontsize',20,...
        'Xtick',yearlength*floor(linspace(0,max(n-1,1),min(n,5))),...
        'XtickLabels',floor(linspace(0,max(n-1,1),min(n,5))),...
        'Yscale',yscl,'Ylim',ylims{j},'Ytick',yticks{j},'Xlim',[0 yearlength*max(n-1,1)])
    grid on;

    if isequal(yscl,'log')
        legend(h2(2),{'predicted data'},'location','sw','interpreter','latex','fontsize',16)
    else
        legend(h2(2),{'predicted data'},'location','nw','interpreter','latex','fontsize',16)
    end

    if j==1
        ylabel(['host $(S,N)$'],'interpreter','latex')
    else
        ylabel(['pathogen $(P,Z)$'],'interpreter','latex')
    end
    % end
    xlabel('generation number ($n$)','interpreter','latex')
    % if j~=1
    %     hold on
    %     h3 = plot(yearlength*(n_err_tol*[1 1]-1),[yl(1) yl(2)],'k--','linewidth',3);
    %     hold off
    %     legend([h1;h2;h3],{'${\bf X_test}$','${\bf Y_test}$','$n_{0.2}(\hat{\bf w})$'},'location','ne','interpreter','latex','fontsize',16)
    % end
    saveas(gcf,['~/Desktop/hybrid_pred',num2str(j),'_snr',num2str(snr_Y),'_numI',num2str(num_train_inds),'_subt,',num2str(subsamp_t),'_rng',num2str(seed1),'.png'])

end

%%

[~,a]=sort(cellfun(@(Y_test)max(Y_test(:,j)),Y_train),'descend');
ks = a(randperm(length(a),3));
ks = (1)';
for j=1:2
figure(j+4);clf
hold on
for k=ks'
    i = train_inds(X_in(k));
    h0=plot(t_epi_test{i},Ycell_test{i}(:,j),'markersize',7,'linewidth',4)
    plot(train_time{i==train_inds(X_in)},...
        Y_train{i==train_inds(X_in)}(:,j)*nY(j),'ro-','markersize',4,'linewidth',3)
    plot(0,X_train(i==train_inds,j)*nX(j),xO_cl,'markersize',12,'linewidth',5);
    % plot(yearlength,X_train(find(i==train_inds)+1,j)*nX(j),xO_cl,'markersize',7,'linewidth',5);
    xlim([0 56])
    grid on;
    set(gca,'ticklabelinterpreter','latex','fontsize',16,'Xtick',[0:7:56],'Yscale','log')
        legend({'true data','observed data'},'location','best','interpreter','latex','fontsize',16)
     set(h0,'color',y_cl)
end
saveas(gcf,['~/Desktop/hybrid_sol_data_zoom',num2str(j),'_',num2str(snr_Y),'_rng',num2str(seed1),'.png'])
hold off
end

