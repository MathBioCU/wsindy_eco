if isequal(varget,'X')
    w = W_X;
    if all([exist('W_IC_true','var'),exist('W_Y_true','var'),exist('W_X_true','var')])
        w_true = W_X_compare;
    else
        w_true = [];
    end
    C = CovW_X;
elseif isequal(varget,'Y')
    w = W_Y;
    if all([exist('W_IC_true','var'),exist('W_Y_true','var'),exist('W_X_true','var')])
        w_true = W_Y_compare;
    else
        w_true = [];
    end
    C = CovW_Y;
elseif isequal(varget,'IC')
    w = W_IC;
    if all([exist('W_IC_true','var'),exist('W_Y_true','var'),exist('W_X_true','var')])
        w_true = W_IC_compare;
    else
        w_true = [];
    end
    C = CovW_IC;
end    
[w_hat,ss] = wendy_param(w);
[w_hat_true,ss_true] = wendy_param(w_true);
xflip = [1:length(w_hat) length(w_hat):-1:1];
c = 0.05; % <(100)c chance of not containing true val
stdW = max(sqrt(diag(C)),eps);
conf_int = arrayfun(@(x)norminv(1 - c/2,0,x),stdW);
h1=plot(1:length(w_hat),w_hat,'ro');hold on
extra_true = [];
ms = cellfun(@(s)length(s),ss);
for j=1:length(w_hat)
    gg = 0.4;
    fill([j-gg j+gg j+gg j-gg],[w_hat(j)-conf_int(j) w_hat(j)-conf_int(j)  w_hat(j)+conf_int(j) w_hat(j)+conf_int(j)],...
        'w','linewidth',1.5);
    line([j-gg j+gg],[w_hat(j) w_hat(j)],'color','r','linewidth',3)
    if ~isempty(ss_true)

        ind = find(j<=cumsum(ms),1);
        [~,jj] = ismember(ss{ind}(j-sum(ms(1:ind-1))),ss_true{ind});

        if jj>0
            plot(j,w_hat_true(sum(length(cell2mat(ss_true(1:ind-1))))+ jj),'*','color','green','linewidth',3,'markersize',15)
        end
    end
end

L = length(w_hat);
for k=1:length(ss)
    if ~isempty(ss_true)
        a1 = setdiff(ss_true{k},ss{k});
    else
        a1 = [];
    end
    for j=L+1:L+length(a1)
        line([j-gg j+gg],[0 0],'color','r','linewidth',3)
        if ~isempty(w_true)
            plot(j,w_true{k}(a1(j-L)),'*','color','green','linewidth',3,'markersize',15)
        end
    end
    L = L + length(a1);
end
% a2 = setdiff(ss_true{2},ss{2});
% for j=length(w_hat)+length(a1)+1:length(w_hat)+length(a1)+length(a2)
%     line([j-gg j+gg],[0 0],'color','r','linewidth',3)
%     plot(j,w_true{2}(a2(j-length(w_hat)-length(a1))),'*','color','green','linewidth',3,'markersize',15)
% end

hold off
hh = get(gca,'children');
if ~isempty(w_true)
    legend(hh([3 2 1]),{[num2str((1-c)*100),'% CI'],'learned val.','true val.'},'fontsize',20,'location','bestoutside','interpreter','latex')
else
    legend(hh([2 1]),{[num2str((1-c)*100),'% CI'],'learned val.'},'fontsize',20,'location','bestoutside','interpreter','latex')
end
set(gca,'Xtick',1:length(w_hat))
xlabel(['$\widehat{\bf w}^{',varget,'}$'],'interpreter','latex')
xlim([0 length(w_hat)+1])
grid on

set(gca,'Xticklabels',arrayfun(@(i)['$w_',num2str(i),'$'],1:length(w_hat),'Un',0))
% xlabel('Host density ($N_n$)','interpreter','latex')
% ylabel('Generation ($n$)','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','fontsize',20)
% saveas(gcf,['~/Desktop/conf_int_',varget,'_',num2str(snr_Y),'.png'])
