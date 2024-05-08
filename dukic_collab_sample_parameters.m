%% get samples
n= 200;
X_pred_cloud = cell(n,1);
x0 = X_test(1,:);

parfor j=1:n
    rng('shuffle')
    disp(j)
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

%% quick view 

figure(1);clf
hold on
cellfun(@(x)semilogy(x(:,1),'b'),X_pred_cloud)
plot(X_test(:,1),'linewidth',3,'color','red')
plot(X_pred(:,1),'linewidth',3,'color','green')
set(gca,'Yscale','log')
legend('X_1')
hold off
saveas(gcf,['~/Desktop/X1_cloud_',num2str(snr_Y),'.png'])

figure(2);clf
hold on
cellfun(@(x)semilogy(x(:,2),'b'),X_pred_cloud)
plot(X_test(:,2),'linewidth',3,'color','red')
plot(X_pred(:,2),'linewidth',3,'color','green')
set(gca,'Yscale','log')
legend('X_2')
hold off

%% get inter-peak distributions

close all;
clz = {[1 0.6 0.3],[0 0.6 1],[0.6 0.6 0.6]};
legs = {{'IPD distrib','True IPD','Mean IPD'}, {'Amp distrib','True Amp','Mean Amp'} };
for ind = 1:2
IPD = [];
Amps = [];
max_gen_cap = 80;
[A_true,IPD_true] = findpeaks(X(1:max_gen_cap,ind));
IPD_true = mean(diff(IPD_true));
A_true = mean(A_true);
for j=1:n
    if length(X_pred_cloud{j}(:,ind))>=3
        [A,I] = findpeaks(X_pred_cloud{j}(1:min(max_gen_cap,end),ind));
        IPD = [IPD;diff(I)];
        Amps = [Amps;A];
    end
end
if ind==1
figure(1+2*(ind-1));
[h,xi] = ksdensity(IPD);
fill(xi,h,clz{3},'linewidth',2,'FaceAlpha',0.4) ;
hold on
plot(IPD_true*[1 1],[0 max(h)],'linewidth',2)
plot(mean(IPD)*[1 1],[0 max(h)],'linewidth',2)
legend(legs{1})
hold off
end

figure(2+2*(ind-1));
[h,xi] = ksdensity(Amps);
fill(xi,h,clz{ind},'linewidth',2,'FaceAlpha',0.4) ;
hold on
plot(A_true*[1 1],[0 max(h)],'linewidth',2)
plot(mean(Amps)*[1 1],[0 max(h)],'linewidth',2)
legend(legs{2})
hold off
end

%% cloud plot
dr = '~/Dropbox/Boulder/research/data/dukic collab/';
load([dr,'UQ_plots_correct_model_snry05.mat'])
ind = 1;
X_ext = cellfun(@(X)[X;NaN*ones(81-size(X,1),nstates_X)],X_pred_cloud,'Un',0);
clf
xbins = linspace(-7,0,50);
Nmin = 20;
Nmax = 60;%size(X_pred,1);
upper_thresh =100; 
figure(1);clf
hold on
cellfun(@(x)plot(log10(x(1:size(x,1),ind)),0:size(x,1)-1,'linewidth',0.5,'color',[0.8 0.8 0.8]),X_pred_cloud)
gap = 2;
for j=Nmax:-gap:Nmin
    dat = cellfun(@(X)X(j,ind),X_ext);
    dat = dat(~isnan(dat));
    dat = dat(dat>=0);
    dat = log10(dat);
    % dat = dat(dat<upper_thresh);
    % xi = linspace(nanmean(dat)-4*nanstd(dat),nanmean(dat)+4*nanstd(dat),50);
    [h,xi] = ksdensity(dat(~isnan(dat)));
    % [h,e] = histcounts(log10(dat),xbins,'normalization','pdf');
    fill(xi,h/9*gap*9+j-1,[1 0.6 0.3],'linewidth',2) ;
    hold on
end
hh=plot(log10(X_test(:,ind)),0:80,'go-',log10(X_pred(1:Nmax+1,ind)),0:Nmax,'b--','linewidth',4,'markersize',12);
% plot(log10(X(:,1)),0:80,'ro-','linewidth',4,'markersize',12)
legend(hh,{'true model output','learned model output'},'location','sw','interpreter','latex');
ylim([Nmin-1 Nmax])
% grid on
view([90 -90])


set(gca,'Xtick',-16:2:16,'Xticklabels',num2str(10.^(-16:2:16)','%0.0e'),'Xlim',[-14 6])
if ind == 1
    xlabel('Host density ($N_n$)','interpreter','latex')
else
    xlabel('Path. density ($Z_n$)','interpreter','latex')
end
ylabel('Generation ($n$)','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','fontsize',18)
set(gcf,'position',[1970         232        1249         689])
% grid on

% saveas(gcf,['~/Desktop/cloud_density_snrY01.png'])
%% view conf int

varget = 'Y';
view_conf_int;
