dr = '~/Desktop/mits_study/';%
loadvars = {'results_cell','ntrain_inds','rngs','sim_cell'};
peaks_ = true;
subt = 2;
ttf = [0.75];
mits = 0;
res_ind = {};
mns = []; meds = [];
for kk=1:6
    snr_Y = 0.1*2^(kk-6);
    sind = 11;%[1 3 4 6 7 9 11]
    load([dr,'sweep_snrY_',num2str(snr_Y),'_ttf_',num2str(ttf),'_subt_',num2str(subt),'_mits_',num2str(mits),'_peaks.mat'],loadvars{:})
    err_tol = 0.5;
    n_err_tols = cellfun(@(s) get_n_err_tol(s{2}{1},s{1}{1},err_tol) , sim_cell(end,:));
    res_ind = [res_ind,{n_err_tols(:)}];
    mns(kk) = mean(n_err_tols');
    meds(kk) = median(n_err_tols');
end
x = 0.1*2.^(-5:0);

violin(res_ind,'kernel','box',...
'bw',range(cell2mat(res_ind(:)))/size(results_cell,2)^(1/2),'support',[min(cell2mat(res_ind(:)))-eps*range(cell2mat(res_ind(:))) max(cell2mat(res_ind(:)))+eps*range(cell2mat(res_ind(:)))],...
'facecolor',[0 0.5 1],'plotlegend',0,'linewidth',1.7,'x',log10(x));
ylabel(['Generations predicted ($n_{',num2str(err_tol),'}(\hat{\bf w})$)'],'interpreter','latex')
xlabel('Noise level ($\sigma_{NR}$)','interpreter','latex')
set(gca,'Xticklabels',x,'Xtick',log10(x),'Xlim',[min(log10(x))-mean(diff(log10(x)))/2 max(log10(x))+mean(diff(log10(x)))/2],...
    'ticklabelinterpreter','latex','fontsize',18)
set(gca,'Ylim',[0 80])
grid on
saveas(gcf,['~/Desktop/n_err_tol_mits_',num2str(mits),'.png'])
