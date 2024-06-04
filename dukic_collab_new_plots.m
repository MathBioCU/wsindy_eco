%% Calculate AIC
dr = '~/Dropbox/Boulder/research/data/dukic collab/';
% load([dr,'UQ_plots_uncorrected_red_model.mat'],'WS_IC','WS_Xeq','WS_Yeq')
load([dr,'UQ_plots_corrected_red_model.mat'],'WS_IC','WS_Xeq','WS_Yeq')
r = 0;num_param=0;diag_reg = 10^-1;
[G_IC,b_IC]=WS_IC.apply_cov(WS_IC.G{1},WS_IC.b{1},diag_reg);
C = WS_IC.cov;
r = r + norm(G_IC*WS_IC.weights-b_IC)^2; %+ log(det((1-diag_reg)*C + diag_reg*speye(size(C,1))));
num_param = num_param + length(find(WS_IC.weights));

[G_Y,b_Y]=WS_Yeq.apply_cov(WS_Yeq.G{1},WS_Yeq.b{1},diag_reg);
C = WS_Yeq.cov;
r = r + norm(G_Y*WS_Yeq.weights-b_Y)^2; %+ log(det((1-diag_reg)*C + diag_reg*speye(size(C,1))));
num_param = num_param + length(find(WS_Yeq.weights));

[G_X,b_X]=WS_Xeq.apply_cov(WS_Xeq.G{1},WS_Xeq.b{1},diag_reg);
C = WS_Xeq.cov;
r = r + norm(G_X*WS_Xeq.weights-b_X)^2; %+ log(det((1-diag_reg)*C + diag_reg*speye(size(C,1))));
num_param = num_param + length(find(WS_Xeq.weights));
r
2*num_param
AIC = 2*num_param+r

%%
dr = '~/Dropbox/Boulder/research/data/dukic collab/';
% load([dr,'UQ_plots_snrY01.mat'])
load([dr,'UQ_plots_correct_model_snry05.mat'])c
% load([dr,'UQ_plots_uncorrected_red_model.mat'])
% load([dr,'UQ_plots_corrected_red_model.mat'])

%% get inter-peak distributions
nstates =2;
PD = cell(n,2);
AM = cell(n,2);
len = zeros(n,2);
for ind = 1:2
    max_gen_cap = 80;
    % [A_true,IPD_true] = findpeaks(X(1:max_gen_cap,ind));
    % IPD_true = mean(diff(IPD_true));
    % A_true = mean(A_true);
    for j=1:n
        if length(X_pred_cloud{j}(:,ind))>=3
            [AM{j,ind},PD{j,ind}] = findpeaks(X_pred_cloud{j}(1:min(max_gen_cap,end),ind),'MinPeakDistance',4);
            len(j,ind) = length(AM{j,ind});
        end
    end
end

mus_P = []; mus_A = [];
for ind = 1:nstates
for peak_num = 1:max(len(:))
    dat = cellfun(@(t)t(peak_num),PD(len(:,ind)>=peak_num,ind));
    % mus_P(peak_num,ind)= (mean(dat)-1)*yearlength;
    stds_P(peak_num,ind)= std(dat)*yearlength;
    mus_P(peak_num,ind)= (10.^mean(log10(dat))-1)*yearlength;
    % stds_P(peak_num,ind)= 10.^std(log10(dat))*yearlength;

    dat = cellfun(@(t)t(peak_num),AM(len(:,ind)>=peak_num,ind));
    % mus_A(peak_num,ind)= mean(dat);
    stds_A(peak_num,ind)= std(dat);
    mus_A(peak_num,ind)= 10.^mean(log10(dat));
    % stds_A(peak_num,ind)= 10.^std(log10(dat));
end
end

%% examine specific peak stats

peak_num = 5;
ind = 1;
dat_A = cellfun(@(t)t(peak_num),AM(len(:,ind)>=peak_num,ind));
dat_P = cellfun(@(t)t(peak_num),PD(len(:,ind)>=peak_num,ind));
scatter(dat_P,dat_A)

%% plot learned trajectory, CI

n = size(X_pred,1);
x_cl = 'ko'; % colors for true data
y_cl = 1.3*[0.4660 0.6740 0.1880]; % colors for true data
xL_cl = 'ko'; % colors for learned data
yL_cl = 'b-';% colors for learned data
xO_cl = 'ko';% colors for observed data
yO_cl = 'r'; % colors for observed data
yobs_cl = 'r-';
yscl = 'log';

err_tol = 0.5;
n_err_tol = get_n_err_tol(X_pred,X_test,err_tol);

for j=1:nstates
    figure(j);clf; set(gcf,'position',[232   409   990 474])
    h2=plot(tn_pred(1:n),X_pred(1:n,j),xL_cl,...
        t_pred, Y_pred(:,j),yL_cl,yearlength*[n_err_tol]*[1 1],ylims{j},'k--','linewidth',3,'markersize',7);
    hold on
    h0=plot(tn_test(1:n),X_test(1:n,j),x_cl,...
        t_test, Y_test(:,j),'--',yearlength*[n_err_tol]*[1 1],ylims{j},'k--','linewidth',3,'markersize',7);
     set(h0(2),'color',y_cl)
    num_std = 2;
    for k=1:size(mus_A,1)
        % h3=rectangle('position',[mus_P(k,j)-num_std*stds_P(k,j) mus_A(k,j)-num_std*stds_A(k,j) 2*num_std*stds_P(k,j) 2*num_std*stds_A(k,j)],...
        %     'facecolor',[1 0 0 0.5])
        h3=fill([mus_P(k,j)-num_std*stds_P(k,j) mus_P(k,j)+num_std*stds_P(k,j) mus_P(k,j)+num_std*stds_P(k,j) mus_P(k,j)-num_std*stds_P(k,j)],...
            [mus_A(k,j)-num_std*stds_A(k,j) mus_A(k,j)-num_std*stds_A(k,j) mus_A(k,j)+num_std*stds_A(k,j) mus_A(k,j)+num_std*stds_A(k,j)],...
            'r','facealpha',0.25);
    end
    % h3=errorbar(mus_P(:,j),mus_A(:,j),-num_std*stds_A(:,j),num_std*stds_A(:,j),...
        % -num_std*stds_P(:,j),num_std*stds_P(:,j),'ro','linewidth',3,'CapSize',15);
    scatter(mus_P(:,j),mus_A(:,j),50,[1 1 0],'filled')
    % scatter(mus_P(1:end-1,j),mus_A(1:end-1,j),50,[1 1 0],'filled')
    set(gca,'ticklabelinterpreter','latex','fontsize',20,...
        'Xtick',yearlength*floor(linspace(0,max(n-1,1),min(n,5))),...
        'XtickLabels',floor(linspace(0,max(n-1,1),min(n,5))),...
        'Yscale',yscl,'Ylim',ylims{j},'Ytick',yticks{j},'Xlim',[0 yearlength*max(n-1,1)])
    grid on;

    legend([h0(2);h2(2);h3],{'true model output','learned model output','95% CR'},'location','sw','interpreter','latex','fontsize',16)
    if j==1
        ylabel(['host $(S,N)$'],'interpreter','latex')
    else
        ylabel(['pathogen $(P,Z)$'],'interpreter','latex')
    end
    xlabel('generation number ($n$)','interpreter','latex')
    saveas(gcf,['~/Desktop/hybrid_pred',num2str(j),'_snr',num2str(snr_Y),'_numI',num2str(num_train_inds),'_subt,',num2str(subsamp_t),'_rng',num2str(seed1),'.png'])
end

%% plot single gen noisy data

[~,a]=sort(cellfun(@(Y_test)max(Y_test(:,j)),Y_train),'descend');
ks = a(randperm(length(a),3));
ks = (1)';
for j=1:2
figure(j+5);clf
hold on
for k=ks'
    i = train_inds(X_in(k));
    h0=plot(t_epi_test{i},Ycell_test{i}(:,j),'markersize',7,'linewidth',4);
    plot(train_time{i==train_inds(X_in)},...
        Y_train{i==train_inds(X_in)}(:,j)*nY(j),'ro-','markersize',4,'linewidth',3)
    plot(0,X_train(i==train_inds,j)*nX(j),xO_cl,'markersize',12,'linewidth',5);
    % plot(yearlength,X_train(find(i==train_inds)+1,j)*nX(j),xO_cl,'markersize',7,'linewidth',5);
    xlim([0 56])
    grid on;
    xlabel('time (days)','interpreter','latex')
    set(gca,'ticklabelinterpreter','latex','fontsize',16,'Xtick',[0:7:56],'Yscale','log')
        legend({'true model output','observed data'},'location','best','interpreter','latex','fontsize',16)
     set(h0,'color',y_cl)
end
saveas(gcf,['~/Desktop/hybrid_sol_data_zoom',num2str(j),'_',num2str(snr_Y),'_rng',num2str(seed1),'.png'])
hold off
end

%% data

close all
for j=1:nstates_X
    figure(1);clf;
    grid on
    xlabel('generation number ($n$)','interpreter','latex')
    for i=1:size(X_test,1)
        if ismember(i,train_inds)
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
    legend([h2],{'observed data'},'location','sw','interpreter','latex','fontsize',16)
    set(gca,'ticklabelinterpreter','latex','fontsize',20,...
        'Xtick',yearlength*floor(linspace(0,max(n-1,1),min(n,5))),...
        'XtickLabels',floor(linspace(0,max(n-1,1),min(n,5))),...
        'Yscale',yscl,'Ytick',yticks{j},'Xlim',[0 yearlength*max(train_inds)+3])
    saveas(gcf,['~/Desktop/hybrid_data',num2str(j),'_snr',num2str(snr_Y),'_numI',num2str(num_train_inds),'_subt,',num2str(subsamp_t),'_rng',num2str(seed1),'_data.png'])
end

%% view conf int

vargets = {'IC','Y','X'};
tiledlayout(1,3)
for i=1:3
    nexttile
    varget = vargets{i};
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

% abs(w_hat)./stdW
conf_int

h1=plot(1:length(w_hat),w_hat,'ro');hold on
extra_true = [];
for j=1:length(w_hat)
    gg = 0.4;
    fill([j-gg j+gg j+gg j-gg],[w_hat(j)-conf_int(j) w_hat(j)-conf_int(j)  w_hat(j)+conf_int(j) w_hat(j)+conf_int(j)],...
        'w','linewidth',1.5);
    if j<=length(ss{1})
        [~,jj] = ismember(ss{1}(j),ss_true{1});
        if jj>0
            plot(j,w_hat_true(jj),'*','color','green','linewidth',2.3,'markersize',10)
        else
            plot(j,0,'*','color','green','linewidth',3,'markersize',10)
        end
    else
        [~,jj] = ismember(ss{2}(j-length(ss{1})),ss_true{2});
        if jj>0
            plot(j,w_hat_true(length(ss_true{1})+jj),'*','color','green','linewidth',2.3,'markersize',10)
            % legend('true vals')
        else
            plot(j,0,'*','color','green','linewidth',2.3,'markersize',10)
        end        
    end
    line([j-gg j+gg],[w_hat(j) w_hat(j)],'color','r','linewidth',3)
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
set(gca,'Xtick',1:length(w_hat))
xlabel(['$\widehat{\bf w}^{',varget,'}$'],'interpreter','latex')
xlim([0 length(w_hat)+1])
grid on

set(gca,'Xticklabels',arrayfun(@(i)['$w_',num2str(i),'$'],1:length(w_hat),'Un',0))
% xlabel('Host density ($N_n$)','interpreter','latex')
% ylabel('Generation ($n$)','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','fontsize',16)
end
hh = get(gca,'children');
legend(hh([2 1 3]),{'true val.','learned val.',[num2str((1-c)*100),'% CI']},'location','bestoutside','interpreter','latex')
set(gcf,'pos',[900 539 1100 250]);
saveas(gcf,['~/Desktop/conf_int_',num2str(snr_Y),'.png'])


