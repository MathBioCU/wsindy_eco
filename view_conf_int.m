if isequal(varget,'X')
w = W_X;
w_true = W_X_compare;
C = CovW_X;
elseif isequal(varget,'Y')
w = W_Y;
w_true = W_Y_compare;
C = CovW_Y;
elseif isequal(varget,'IC')
w = W_IC;
w_true = W_IC_compare;
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
for j=1:length(w_hat)
    gg = 0.4;
    fill([j-gg j+gg j+gg j-gg],[w_hat(j)-conf_int(j) w_hat(j)-conf_int(j)  w_hat(j)+conf_int(j) w_hat(j)+conf_int(j)],...
        'w','linewidth',1.5);
    line([j-gg j+gg],[w_hat(j) w_hat(j)],'color','r','linewidth',3)
    if j<=length(ss{1})
        [~,jj] = ismember(ss{1}(j),ss_true{1});
        if jj>0
            plot(j,w_hat_true(jj),'*','color','green','linewidth',3,'markersize',15)
        end
    else
        [~,jj] = ismember(ss{2}(j-length(ss{1})),ss_true{2});
        if jj>0
            plot(j,w_hat_true(length(ss_true{1})+jj),'*','color','green','linewidth',3,'markersize',15)
            % legend('true vals')
        end        
    end
end

a1 = setdiff(ss_true{1},ss{1});
for j=length(w_hat)+1:length(w_hat)+length(a1)
    line([j-gg j+gg],[0 0],'color','r','linewidth',3)
    plot(j,w_true{1}(a1(j-length(w_hat))),'*','color','green','linewidth',3,'markersize',15)
end
a2 = setdiff(ss_true{2},ss{2});
for j=length(w_hat)+length(a1)+1:length(w_hat)+length(a1)+length(a2)
    line([j-gg j+gg],[0 0],'color','r','linewidth',3)
    plot(j,w_true{2}(a2(j-length(w_hat)-length(a1))),'*','color','green','linewidth',3,'markersize',15)
end

hold off
hh = get(gca,'children');
legend(hh([3 2 1]),{[num2str((1-c)*100),'% CI'],'learned val.','true val.'},'fontsize',20,'location','bestoutside','interpreter','latex')
set(gca,'Xtick',1:length(w_hat))
xlabel(['$\widehat{\bf w}^{',varget,'}$'],'interpreter','latex')
xlim([0 length(w_hat)+1])
grid on

set(gca,'Xticklabels',arrayfun(@(i)['$w_',num2str(i),'$'],1:length(w_hat),'Un',0))
% xlabel('Host density ($N_n$)','interpreter','latex')
% ylabel('Generation ($n$)','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','fontsize',20)
saveas(gcf,['~/Desktop/conf_int_',varget,'_',num2str(snr_Y),'.png'])
