figure(2);clf
x_cl = 'ko'; % colors for true data
y_cl = 1.3*[0.4660 0.6740 0.1880]; % colors for true data
xL_cl = 'ko'; % colors for learned data
yL_cl = 'b--';% colors for learned data
xO_cl = 'ko';% colors for observed data
yO_cl = 'r'; % colors for observed data
yobs_cl = 'r-';
for j=1:nstates_X
    subplot(1,nstates_X,j)
    h0=plot(tn_test(1:n),X_test_cell{jj}(1:n,j),x_cl,...
        t_test, Y_test(:,j),'-',yearlength*[n_err_tol]*[1 1],ylims{j},'k--','linewidth',3,'markersize',7);
    set(h0(2),'color',y_cl)
    hold on
    h2=plot(tn_pred(1:n),X_pred_cell{jj}(1:n,j),xL_cl,...
        t_pred, Y_pred(:,j),yL_cl,yearlength*[n_err_tol]*[1 1],ylims{j},'k--','linewidth',3,'markersize',7);
    set(gca,'ticklabelinterpreter','latex','fontsize',12,...
        'Xtick',yearlength*floor(linspace(0,max(n-1,1),min(n,5))),...
        'XtickLabels',floor(linspace(0,max(n-1,1),min(n,5))),...
        'Yscale',yscl,'Ylim',ylims{j},'Ytick',yticks{j},'Xlim',[0 yearlength*max(n-1,1)])
    grid on;

    legend([h0(2);h2(2)],{'true model output','learned model output'},'location','sw','interpreter','latex','fontsize',14)
    if j==1
        ylabel(['host $(S,N)$'],'interpreter','latex')
    else
        ylabel(['pathogen $(P,Z)$'],'interpreter','latex')
    end
    xlabel('generation number ($n$)','interpreter','latex')  
    % xlim([0 20*yearlength])
end
%%
figure(3);clf
for j=1:2
    subplot(1,2,j)
    h0=plot(t_test, Y_test(:,3+j),'linewidth',3,'markersize',7);
    set(h0(1),'color',    y_cl)
    hold on
    h2=plot(...
        t_pred, Y_pred(:,3+j),yL_cl,'linewidth',3,'markersize',7);
    set(gca,'ticklabelinterpreter','latex','fontsize',12,...
        'Xtick',yearlength*floor(linspace(0,max(n-1,1),min(n,5))),...
        'XtickLabels',floor(linspace(0,max(n-1,1),min(n,5))),...
        'Yscale',yscl,'Xlim',[0 yearlength*max(n-1,1)])
    grid on;
    
    legend([h0(1);h2(1)],{'true model output','learned model output'},'location','sw','interpreter','latex','fontsize',14)
    ylabel(['susc. ',num2str(j),' $(\nu_',num2str(j),')$'],'interpreter','latex')
    xlabel('generation number ($n$)','interpreter','latex') 
    xlim([0 20*yearlength])
end
%%

j=1;
figure(4);clf
i=1;
subplot(3,ceil(size(Y_ns,1)/3+1/3),i)
X_sub = find(diff(train_inds)==1);
plot(train_time{i},nY(j)*Y_train{X_in==X_sub(i)}(:,j),'r-',Y_ns{i,2},Y_ns{i,1}(:,j),'b-.',...
    'linewidth',3)
set(gca,'fontsize',40,'Xtick', [],'Ytick',[])
legend({'data','learned'},'location','westoutside')
for i=1:size(Y_ns,1)
    subplot(3,ceil(size(Y_ns,1)/3+1/3),i+1)
    if ismember(X_sub(i),X_in)
        plot(train_time{i},nY(j)*Y_train{X_in==X_sub(i)}(:,j),'r-',Y_ns{i,2},Y_ns{i,1}(:,j),'b-.',...
            'linewidth',3)
        set(gca,'fontsize',14)
        if i>=2*ceil(size(Y_ns,1)/3+1/3)
            xlabel('t')
        end
    end
end
sgtitle('Host (S)')