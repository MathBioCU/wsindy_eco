%% sweep
addpath(genpath('wsindy_obj_base'))
snr_Ys = 0.01;
ntrain_inds = 9;
rngs = 1:200;

snr_X = 0; %<<< sweep over 
subsamp_t = 1; %<<< sweep over

noise_alg_X = 'logn'; noise_alg_Y = 'logn'; %<<< fixed
test_length = 40; %<<< fixed
err_tol = 0.2; %<<< fixed
stop_tol = 10; %<<< fixed
toggle_zero_crossing = 1; %<<< fixed

tol_dd_sim = 10^-10; %<<< doesn't affect alg
phifun_Y = optTFcos(3,0);%@(t)exp(-9./(1-t.^2));%<<< user choice
tf_Y_params = {'meth','direct','param',9,'mtmin',5,'subinds',3};%<<< user choice
WENDy_args = {'maxits_wendy',5,...
    'lambdas',10.^linspace(-4,0,50),...
    'ittol',10^-4,'diag_reg',10^-6,'verbose',0};
autowendy = 1; %<<< decision needs to made
tol = 5; %<<< decision needs to made
tol_min = 0.1; %<<< decision needs to made
tol_dd_learn = 10^-8;%<<< decision made

pmax_IC = 4;%<<< decision made 
polys_Y_Yeq = 0:3; %<<< decision needs to made
polys_X_Xeq = 0:2; %<<< decision needs to made
pmax_X_Yeq = 4; %<<< decision made 
pmax_Y_Xeq = 4; %<<< decision made
neg_Y = 0; %<<< decision needs to made
neg_X = 0; %<<< decision needs to made
boolT = @(T)all([min(T,[],2)>=-2 sum(T,2)>=-2 max(T,[],2)<4],2); %<<< decision made

dr = '/home/danielmessenger/Dropbox/Boulder/research/data/dukic collab/';
% load([dr,'FH_feedback.mat']);
load([dr,'Gregs_mod_V=0.5.mat'],'Ycell','X','t_epi','custom_tags_X',...
    'yearlength','custom_tags_Y','linregargs_fun_IC','linregargs_fun_Y',...
    'linregargs_fun_X','nstates_X','nstates_Y','W_IC_true','tags_IC_true',...
    'W_Y_true','tags_X_true','tags_Y_true','W_X_true','tags_Ext_X_true','tags_Ext_Y_true','tn','t','Y');

results_cell = cell(length(ntrain_inds),length(rngs));
for train_time_frac = 0.25 %<<< sweep over
for kk=1:length(snr_Ys)
    snr_Y = snr_Ys(kk);
    for ii=1:length(ntrain_inds)
        for jj=1:length(rngs)
             disp([ii jj])
             rng(rngs(jj)); rng_seed = rng().Seed; rng(rng_seed);
    
            num_train_inds = ntrain_inds(ii);      
            %%% get data
            [Y_train,X_train,train_inds,train_time,nstates_X,nstates_Y,X_in,sigma_X,sigma_Y,nX,nY] = ...
                format_data(Ycell,X,t_epi,subsamp_t,train_time_frac,num_train_inds,test_length,snr_X,snr_Y,noise_alg_X,noise_alg_Y);
            toggle_sim = size(X,1);
    
            tic,
            %%% run alg
            [rhs_IC,W_IC,rhs_Y,W_Y,rhs_X,W_X,lib_Y_IC,lib_X_IC,lib_Y_Yeq,lib_X_Yeq,lib_Y_Xeq,lib_X_Xeq]= ...
                wsindy_eco_fcn(toggle_zero_crossing,stop_tol,phifun_Y,tf_Y_params,WENDy_args,autowendy,tol,tol_min,tol_dd_learn,pmax_IC,polys_Y_Yeq,polys_X_Xeq,pmax_X_Yeq,pmax_Y_Xeq,neg_Y,neg_X,boolT,...
                Y_train,X_train,train_inds,train_time,t_epi,yearlength,nstates_X,nstates_Y,X_in,nX,nY,...
                custom_tags_X,custom_tags_Y,linregargs_fun_IC,linregargs_fun_Y,linregargs_fun_X);
            RT = toc;
    
            %%% process results
            W_IC_compare = inject_coeff_param(W_IC_true,zeros(1,nstates_Y),tags_IC_true,cell2mat(lib_Y_IC.tags'),cell2mat(lib_X_IC.tags'));
            errs_2_IC = norm(reshape([W_IC{:}]-[W_IC_compare{:}],[],1))/norm(reshape([W_IC_compare{:}],[],1));
            errs_inf_IC = abs(reshape([W_IC{:}]-[W_IC_compare{:}],[],1))./abs(reshape([W_IC_compare{:}],[],1));
            errs_inf_IC = max(errs_inf_IC(errs_inf_IC<inf),[],'all','omitnan');
            if isempty(errs_inf_IC)
                errs_inf_IC = NaN;
            end
            tpr_IC = tpscore(reshape([W_IC{:}],[],1),reshape([W_IC_compare{:}],[],1));
            
            W_Y_compare = inject_coeff_param(W_Y_true,tags_Y_true,tags_X_true, cell2mat(lib_Y_Yeq.tags'),cell2mat(lib_X_Yeq.tags'));
            errs_2_Y = norm(reshape([W_Y{:}]-[W_Y_compare{:}],[],1))/norm(reshape([W_Y_compare{:}],[],1));
            errs_inf_Y = abs(reshape([W_Y{:}]-[W_Y_compare{:}],[],1))./abs(reshape([W_Y_compare{:}],[],1));
            errs_inf_Y = max(errs_inf_Y(errs_inf_Y<inf),[],'all','omitnan');
            if isempty(errs_inf_Y)
                errs_inf_Y = NaN;
            end
            tpr_Y = tpscore(reshape([W_Y{:}],[],1),reshape([W_Y_compare{:}],[],1));
            
            W_X_compare = inject_coeff_param(W_X_true,tags_Ext_X_true,tags_Ext_Y_true,cell2mat(lib_X_Xeq.tags'),cell2mat(lib_Y_Xeq.tags'));
            errs_2_X = norm(reshape([W_X{:}]-[W_X_compare{:}],[],1))/norm(reshape([W_X_compare{:}],[],1));
            errs_inf_X = abs(reshape([W_X{:}]-[W_X_compare{:}],[],1))./abs(reshape([W_X_compare{:}],[],1));
            errs_inf_X = max(errs_inf_X(errs_inf_X<inf),[],'all','omitnan');
            if isempty(errs_inf_X)
                errs_inf_X = NaN;
            end
            tpr_X = tpscore(reshape([W_X{:}],[],1),reshape([W_X_compare{:}],[],1));
    
            num_gen = floor(size(X,1));
            Ycell_pred = {};
            X_pred = X(1,:);
            options_ode_sim = odeset('RelTol',tol_dd_learn,'AbsTol',tol_dd_sim*ones(1,nstates_Y),'Events',@(T,Y)myEvent(T,Y,5*max(max(cell2mat(Y_train))),toggle_zero_crossing));
            tn_pred = [];
            n=1;
            check=1;
            while and(n<toggle_sim,check)
                rhs_learned = @(y)rhs_Y(y,X_pred(n,:));
                Y0 = rhs_IC(X_pred(n,:));
                t_train = t_epi{n};
                [t_n,Y_n,TE]=ode15s(@(t,x)rhs_learned(x),t_train,Y0,options_ode_sim);
                Ycell_pred = [Ycell_pred;{Y_n}];
                tn_pred = [tn_pred;t_n+(n-1)*yearlength];
                X_pred(n+1,:) = rhs_X(X_pred(n,:),Y_n(end,:));
                if toggle_zero_crossing==1
                    check = all([isempty(TE) X_pred(n+1,:)>=0]);
                else
                    check = isempty(TE);
                end
                n=n+1;
            end
            cumerr = arrayfun(@(i)norm(vecnorm(X_pred(1:i,:)-X(1:i,:),2,2))/norm(vecnorm(X(1:i,:),2,2)),(1:n)');
            n_err_tol = find(cumerr>err_tol,1);
            if isempty(n_err_tol)
                n_err_tol = n-1;
            else
                n_err_tol = n_err_tol-1;
            end
            results_cell{ii,jj} = [errs_2_IC,errs_inf_IC,tpr_IC,...
                errs_2_Y,errs_inf_Y,tpr_Y,...
                errs_2_X,errs_inf_X,tpr_X,RT,n_err_tol];
            % sim_cell{ii,jj} = {X_pred,Ycell_pred};
    
        end
    end
    save([dr,'sweep_snrY_',num2str(snr_Y),'_ttf_',num2str(train_time_frac),'_',date,'.mat'])
end
end
%% view
ylabs = {'Coefficient error ($E_2^{IC}$)',[],'True Positivity Ratio (TPR$^{IC})$','$E_2^{Y}$',[],'TPR$^{Y}$','Coefficient error  ($E_2^{X}$)',[],'True Positivity Ratio (TPR$^{X})$','Walltime(sec)','Generations predicted ($n_{0.2}(\hat{\bf w})$)'};
dr = '~/Dropbox/Boulder/research/data/dukic collab/';
loadvars = {'results_cell','snr_Y','ntrain_inds','rngs'};
subx = 2:9;
for subt = 1
for ttf = [0.5]
for kk = [0.005]
for sind = [7 9 11] %[1 3 4 6 7 9]

disp([sind kk])
% load([dr,'sweep_snrY_',num2str(kk),'_ttf_',num2str(ttf),'_subt_',num2str(subt),'.mat'],loadvars{:}) %ttf = 0.5; ntinds = 6:30
% load([dr,'sweep_snrY_',num2str(kk),'_ttf_',num2str(ttf),'_06-Apr-2024.mat'],loadvars{:})
% load([dr,'sweep_snrY_',num2str(kk),'_ttf_',num2str(ttf),'_07-Apr-2024.mat'],loadvars{:})
% load([dr,'sweep_snrY_',num2str(kk),'_ttf_',num2str(ttf),'_subt_',num2str(subt),'_10-Apr-2024.mat'],loadvars{:})
load([dr,'sweep_snrY_',num2str(kk),'_ttf_',num2str(ttf),'_subt_',num2str(subt),'.mat'],loadvars{:})

runs = length(rngs);
filter_fun = @(r)all([r(3)<=1 r(6)<=1 r(9)<=1]);

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

h=boxplot(res_ind,g,'whisker',1.7);
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

function [value, isterminal, direction] = myEvent(T, Y, thresh, toggle_zero_crossing)
    if toggle_zero_crossing
        value      = or(norm(Y) >= thresh, any(Y==0));
    else
        value      = norm(Y) >= thresh;
    end
    isterminal = 1;   % Stop the integration
    direction  = 0;
end


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
