%% get samples
n=200;
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
plot(X(:,1),'linewidth',3,'color','red')
plot(X_pred(:,1),'linewidth',3,'color','green')
set(gca,'Yscale','log')
legend('X_1')
hold off
saveas(gcf,['~/Desktop/X1_cloud_',num2str(snr_Y),'.png'])

figure(2);clf
hold on
cellfun(@(x)semilogy(x(:,2),'b'),X_pred_cloud)
plot(X(:,2),'linewidth',3,'color','red')
plot(X_pred(:,2),'linewidth',3,'color','green')
set(gca,'Yscale','log')
legend('X_2')
hold off
saveas(gcf,['~/Desktop/X2_cloud_',num2str(snr_Y),'.png'])

%%
X_ext = cellfun(@(X)[X;NaN*ones(81-size(X,1),nstates_X)],X_pred_cloud,'Un',0);
clf
xbins = linspace(-7,0,50);
Nmin = 6;
Nmax = size(X_pred,1);

figure(1);clf
hold on
cellfun(@(x)plot(log10(x(1:size(x,1),1)),0:size(x,1)-1,'linewidth',0.5,'color',[0.8 0.8 0.8]),X_pred_cloud)
gap = 2;
for j=Nmax:-gap:Nmin
    dat = log10(cellfun(@(X)X(j,1),X_ext));
    xi = linspace(mean(dat)-4*std(dat),mean(dat)+4*std(dat),50);
    [h,xi] = ksdensity(dat);
    if j==Nmin
        ff = max(h);
    end
    % [h,e] = histcounts(log10(dat),xbins,'normalization','pdf');
    fill(xi,h/ff*gap*20+j-1,[1 0.6 0.3],'linewidth',2) 
    hold on
end
% plot(log10(X(:,1)),0:80,'ro-',log10(X_pred(:,1)),0:Nmax-1,'g--','linewidth',4,'markersize',12)
plot(log10(X(:,1)),0:80,'ro-','linewidth',4,'markersize',12)
ylim([Nmin 25])
% grid on
view([90 -90])


set(gca,'Xtick',-6:2:2,'Xticklabels',num2str(10.^(-8:2:1)','%0.0e'),'Xlim',[-8 8])
xlabel('Host density ($N_n$)','interpreter','latex')
ylabel('Generation ($n$)','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','fontsize',18)
% grid on

saveas(gcf,['~/Desktop/cloud_density.png'])
%% view conf int
varget = 'Y';
if varget=='X'
w = W_X;
w_true = W_X_compare;
C = CovW_X;
elseif varget=='Y'
w = W_Y;
w_true = W_Y_compare;
C = CovW_Y;
end    
figure(3);clf
[w_hat,ss] = wendy_param(w);
[w_hat_true,ss_true] = wendy_param(w_true);
xflip = [1:length(w_hat) length(w_hat):-1:1];
c = 0.05; % <(100)c chance of not containing true val
stdW = max(sqrt(diag(C)),eps);
conf_int = arrayfun(@(x)norminv(1 - c/2,0,x),stdW);
h1=plot(1:length(w_hat),w_hat,'ro');hold on
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

hold off
hh = get(gca,'children');
legend(hh([3 2 1]),{[num2str((1-c)*100),'% CI'],'learned val.','true val.'},'fontsize',20,'location','bestoutside','interpreter','latex')
set(gca,'Xtick',1:length(w_hat))
xlabel(['$\widehat{\bf w}^',varget,'$'],'interpreter','latex')
xlim([0 length(w_hat)+1])
grid on

set(gca,'Xticklabels',arrayfun(@(i)['$w_',num2str(i),'$'],1:length(w_hat),'Un',0))
% xlabel('Host density ($N_n$)','interpreter','latex')
% ylabel('Generation ($n$)','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','fontsize',24)
saveas(gcf,['~/Desktop/conf_int_',varget,'_',num2str(snr_Y),'.png'])


function [W_p,ss] = wendy_param(W)
    W_p = cell2mat(cellfun(@(w)w(:),W,'Un',0));
    W_p = W_p(W_p~=0);
    ss = cellfun(@(w)find(w),W,'Un',0);
end

function W_s = wendy_hybrid_sample(W,C)
    [W_p,ss] = wendy_param(W);
    W_samp = mvnrnd(W_p,C);
    W_s = W;
    ind = 0;
    for j=1:length(W)
        W_s{j}(ss{j}) = W_samp(ind+1:ind+length(ss{j}));
        ind = ind + length(ss{j});
    end
end
