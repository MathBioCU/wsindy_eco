%% view
addpath(genpath('../utils'))
addpath(genpath('../wsindy_obj_base'))

% dr = '~/Dropbox/Boulder/research/data/dukic collab/';
dr = '~/Desktop/';
peaks_ = true;
subx = 1;
subt = 3;
ttf = [0.75];
nz_levs = [0.002 0.01 0.02 0.05];
err_tol = 0.5;
num_false = 0;
sub_sind = {0,0,0,0,0,0,0};
sinds = [1 3 4 6 7 9 11];
toggle_filter=0;

for iii=0
sub_sind = {iii};
for sinds = [11]

    if iii==0
        sim_inds = 1;
    elseif iii==1
        sim_inds = 2;
    elseif iii==2
        sim_inds = 3:7;
    end

if toggle_filter==0
    filter_fun = @(r)all([r{3}(1:end-1)<=1;r{6}(5)<=1;r{9}(3)<=1]);
elseif toggle_filter==1
    filter_fun = @(r)all([r{3}(1:end-1)<=1;r{6}(4)<=1;r{9}(3)==1]);
end    
loadvars = {'results_cell','snr_Y','ntrain_inds','rngs','sim_cell','num_gen','num_sim'};
ylabs = {'Coefficient error ($E_2^{IC}$)',[],'True Positivity Ratio (TPR$^{IC})$',...
    '$E_2$',[],'TPR',...
    '$E_2$',[],'TPR','Walltime(sec)',...
    ['Gens predicted ($n_{',num2str(err_tol),'}(\hat{\bf w})$)'],...
    ['Peaks predicted']};

for mm = 1:length(sinds)
    sind = sinds(mm);

    rr = {};
    for kk=nz_levs
        load([dr,'TwoPath_snrY_',num2str(kk),'_ttf_',num2str(ttf),'_subt_',num2str(subt),'_mits_5_peaks.mat'],loadvars{:})
        runs = length(rngs);
        for j=1:length(results_cell)
            results_cell{j}{1}=results_cell{j}{1}(1:end-1);
            results_cell{j}{3}=results_cell{j}{3}(1:end-1);
        end
        if sind==11
            n_err_tols = cell2mat(cellfun(@(s)...
                cellfun(@(w,v)get_n_err_tol(w,v,err_tol),s{2}(sim_inds),s{1}(sim_inds)),...
                sim_cell,'un',0));
        elseif sind==12
            n_err_tols = cell2mat(cellfun(@(s)...
                cellfun(@(w,v)peak_diff(w,v,2),s{2}(sim_inds),s{1}(sim_inds)),...
                sim_cell,'un',0));
        end
        length(find(cellfun(@(s)cellfun(@(w)any(w(1:min(10,end),end)<=0),s{2}(sim_inds)),sim_cell)))/runs
        % length(find(cell2mat(cellfun(@(s)cellfun(@(w)any(w(1:min(5,end),end)<=0),s{2}(sim_inds)),sim_cell(:),'un',0))))/runs/5

        filter_inds = cellfun(@(r)filter_fun(r),results_cell);
        disp(['fraction kept=',num2str(arrayfun(@(i)length(find(filter_inds(i,:)))/size(filter_inds,2),1:size(filter_inds,1))')])
        if sind >=11
            res_ind = n_err_tols;
        else
            if isequal(sub_sind{mm},0)
                res_ind = cellfun(@(r)mean(r{sind}),results_cell);
            else
                res_ind = cellfun(@(r)mean(r{sind}(sub_sind{mm})),results_cell);
            end
        end
        num_r = size(res_ind,1);
        res_ind = res_ind(:,filter_inds);
        res_ind = res_ind(:);    
        if ismember(sind,[1 4 7])
            OL = res_ind(res_ind>100);
            disp(['percent remaining errs > 100'])
            res_ind = res_ind(res_ind<=100);
        end
        rr = [rr;{res_ind}];
    end
% figure(sind);clf
% h=histogram(n_err_tols(:)); hold on
% plot(mean(n_err_tols(:)')*[1;1],max(h.Values)*[0;1],'r-',...
%     median(n_err_tols(:)')*[1;1],max(h.Values)*[0;1],'k-')
% legend({'h',['mean=',num2str(mean(n_err_tols(:)'))],['med=',num2str(median(n_err_tols(:)'))]})
% hold off

   
    if num_r>0
        x = nz_levs*100;
        if ~ismember(sind,[1 4 7])
            if range(res_ind)>0
                bw = range(cell2mat(rr))/size(results_cell,2)^(1/2);
            else
                bw = [];
            end
            
            violin(rr','kernel','box',...
            'bw',bw,'support',...
            [min(cell2mat(rr))-10^-10*range(cell2mat(rr)) max(cell2mat(rr))+10^-10*range(cell2mat(rr))],...
            'facecolor',[0 0.5 1],'plotlegend',0,'linewidth',1.7,'x',1:length(nz_levs))%,'x',1:length(nz_levs));
        else
            violin(cellfun(@(r)log10(r),rr','Un',0),'kernel','box',...
            'bw',range(log10(cell2mat(rr)))/size(results_cell,2)^(1/2),'support',[-7 2+eps],...
            'facecolor',[0 0.5 1],'plotlegend',0,'linewidth',1.7,'x',1:length(nz_levs));
        end
        set(gca,'Xlim',[0.5 length(nz_levs)+0.5])
        if length(x)>1
            % set(gca,'Xtick',1:length(nz_levs))%,'Xlim',[min(x)-mean(diff(x))/2 max(x)+mean(diff(x))/2])
        else
            set(gca,'Xtick',x)
        end
        
        ylabel(ylabs(sind),'interpreter','latex')
        xlabel('Noise level ($100\times\sigma_{NR}\%$)','interpreter','latex')   
        set(gca,'Xticklabels',x,'ticklabelinterpreter','latex','fontsize',18)        
        if ismember(sind,[1 4 7])
            set(gca,'Ytick',-6:2:2,'Yticklabels',num2str(10.^(-6:2:2)','%0.0e'),'Ylim',[-6 2])
        elseif ismember(sind,[3 6 9])
            set(gca,'Ylim',[0 1])
        elseif sind==11
            set(gca,'Ylim',[0 num_gen])
        end
        grid on

        if sind==12
            set(gca,'Ytick',0:8)
        end
        c = 0.05; % <(100)c chance of not containing true val
        mus = cellfun(@(r)mean(r),rr(:));
        stdW = cellfun(@(r)std(r),rr(:));        
        conf_int = arrayfun(@(x)norminv(1 - c/2,0,x),stdW);
        disp([mus stdW])
        if toggle_filter==0
            saveas(gcf,['~/Desktop/multipath_stat',num2str(sind),'_eq',num2str(sub_sind{mm}),'.png'])
        elseif toggle_filter==1
            saveas(gcf,['~/Desktop/multipath_stat',num2str(sind),'_eq',num2str(sub_sind{mm}),'_f.png'])
        end
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
violin(res_ind,'kernel','box','bw',500^(1/4),'facecolor',[0 1 1],'plotlegend',0,'linewidth',1.7);
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

% distributionPlot(res_ind','colormap',copper)
ylim([0 80])
%%
sind = 11;
res_ind = cellfun(@(r)r(sind),results_cell(11,:));
figure(1);clf
histogram(res_ind,20)


%%

dr = '~/Dropbox/Boulder/research/data/dukic collab/';
kk = 0.01; ttf = 0.75; subt = 2;
loadvars = {'results_cell','snr_Y','ntrain_inds','rngs','sim_cell','coeffs_cell','libs_cell'};
load([dr,'TwoPath_snrY_',num2str(kk),'_ttf_',num2str(ttf),'_subt_',num2str(subt),'_mits_5_peaks.mat'],loadvars{:})

subx = 2; sind = 9;
res_ind = cellfun(@(r)r(sind),results_cell(subx,:));
inds = find(res_ind<1);
inds_1 = find(res_ind==1);

coeff_ref = coeffs_cell{subx,inds_1(1)}{end};

coeffs = cellfun(@(c)c{end},coeffs_cell(subx,inds),'uni',0);
libs = cellfun(@(c)c(end-1:end),libs_cell(subx,inds),'uni',0);
for coord = 1:2
    coeffs_coord = cellfun(@(c)c{coord},coeffs,'un',0);
    coeff_ref{coord}
    mean(cat(3,coeffs_coord{:}),3)
end
libs{1}{:}
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