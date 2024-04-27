%% view
ylabs = {'Coefficient error ($E_2^{IC}$)',[],'True Positivity Ratio (TPR$^{IC})$','$E_2^{Y}$',[],'TPR$^{Y}$','Coefficient error  ($E_2^{X}$)',[],'True Positivity Ratio (TPR$^{X})$','Walltime(sec)','Generations predicted ($n_{0.2}(\hat{\bf w})$)'};
dr = '~/Dropbox/Boulder/research/data/dukic collab/';
loadvars = {'results_cell','snr_Y','ntrain_inds','rngs'};
subx = ':';
for subt = 1
for ttf = [0.75]
for kk = [0.03]
for sind = [7 9 11] %[1 3 4 6 7 9]

disp([sind kk])
% load([dr,'sweep_snrY_',num2str(kk),'_ttf_',num2str(ttf),'_subt_',num2str(subt),'.mat'],loadvars{:}) %ttf = 0.5; ntinds = 6:30
% load([dr,'sweep_snrY_',num2str(kk),'_ttf_',num2str(ttf),'_06-Apr-2024.mat'],loadvars{:})
% load([dr,'sweep_snrY_',num2str(kk),'_ttf_',num2str(ttf),'_07-Apr-2024.mat'],loadvars{:})
% load([dr,'sweep_snrY_',num2str(kk),'_ttf_',num2str(ttf),'_subt_',num2str(subt),'_10-Apr-2024.mat'],loadvars{:})
load([dr,'sweep_snrY_',num2str(kk),'_ttf_',num2str(ttf),'_subt_',num2str(subt),'.mat'],loadvars{:})

runs = length(rngs);
filter_fun = @(r)all([r(3)<=1 r(6)<=1 r(9)<=1]);

if isequal(subx,':')
    subx = 1:size(results_cell,1);
end

filter_inds = cellfun(@(r)filter_fun(r),results_cell(subx,:));
% arrayfun(@(i)length(find(filter_inds(i,:))),1:size(filter_inds,1))'
res_ind = cellfun(@(r)r(sind),results_cell(subx,:));
if sind==11
    arrayfun(@(i) length(find(res_ind(i,:)>=79)),1:size(res_ind,1))
end
res_ind = arrayfun(@(i)res_ind(i,filter_inds(i,:)),1:length(subx),'uni',0);

if ismember(sind,[1 4 7])
    % res_ind = cellfun(@(r)max(r,10^-6),res_ind,'uni',0);
    OL = cellfun(@(r)r(r>100),res_ind,'uni',0);
    cellfun(@(o)length(o)/runs,OL)
    res_ind = cellfun(@(r)r(r<=100),res_ind,'uni',0);
end
g = arrayfun(@(i) zeros(length(res_ind{i}),1)+i-1,(1:length(res_ind))','uni',0);
g = cell2mat(g);
res_ind = cell2mat(res_ind)';
x = ntrain_inds(subx);
figure(sind);

h=boxplot(res_ind,g,'whisker',1.5);
set(h,'LineWidth',2)
set(h,'MarkerSize',10)
% replaceboxes

% ylim=[0 size(X,1)-1];
% ylabel('$n_{0.2}(\hat{\bf w})$','interpreter','latex')
ylabel(ylabs(sind),'interpreter','latex')
xlabel('Generations observed ($|\mathcal{I}|$)','interpreter','latex')
% legend(h,['$\sigma_{NR}=',num2str(snr_Y),'$'],'interpreter','latex','location','best')
% set(gca,'Xtick',x,'ticklabelinterpreter','latex','fontsize',20,'xlim',[min(x) max(x)])
set(gca,'Xticklabels',x,'ticklabelinterpreter','latex','fontsize',18)
if ismember(sind,[1 4 7])
    set(gca,'Yscale','log','Ytick',10.^(-6:2:2),'Ylim',[10^-6 10^2])
elseif sind==11
    set(gca,'Ylim',[0 80])
end
grid on
saveas(gcf,['~/Desktop/hybrid_snrY_',num2str(kk),'_ttf_',num2str(ttf),'_stat',num2str(sind),'_tpr1.png'])

end
end
end
end

%%
ylabs = {'$E_2^{IC}$',[],'TPR$^{IC}$','$E_2^{Y}$',[],'TPR$^{Y}$','$E_2^{X}$',[],'TPR$^{X}$','Walltime(sec)','$n_{0.2}(\hat{\bf w})$'};
sind = 11;%[1 3 4 6 7 9]
subx = 2:9;
res_ind = cellfun(@(r)r(sind),results_cell(subx,:));
figure(1);clf
% histogram(res_ind(end,:),20)
violin(res_ind')
% distributionPlot(res_ind','colormap',copper)
ylim([0 80])
%%
sind = 11;
res_ind = cellfun(@(r)r(sind),results_cell(11,:));
figure(1);clf
histogram(res_ind,20)


% Q = quantile(res_ind,3,2);
% y = Q(:,2);
% y_l = Q(:,1);
% y_u = Q(:,3);
% y_f = [y_l;flipud(y_u)];
% x_f = [x';flipud(x')];
% fill(x_f,y_f,[0 1 1])
% hold on
% h=plot(x,y,'k-o',x,y_u,'r-',...
%     x,y_l,'r-','linewidth',3);
% hold off
% grid on