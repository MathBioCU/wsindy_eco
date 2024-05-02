addpath(genpath('wsindy_obj_base'))
ntrain_inds = [12 16 20];
rngs = 1:500;

snr_X = 0; % noise level for X
noise_alg_X = 'logn'; % noise distribution for X
noise_alg_Y = 'logn'; % noise distribution for Y

test_length = 40; % number of generations to test over
err_tol = 0.2; % tol for n_tol= number of generations for which cumulative rel err < tol
stop_tol = 100; % halt simulations if values exceed max observed by this multitude
toggle_zero_crossing = 1; % halt simulations that are non-positive

toggle_sim = 1; % toggle perform diagnostic forward simulation
num_sim = 0; % number of out-of-sample testing simulations
oos_std = 0.2; % std of out-of-sample ICs, uniformly randomly sampled around training IC
tol_dd_sim = 10^-10; % ODE tolerance (abs,rel) for diagnostic sim

phifun_Y = @(t)(1-t.^2).^9; % test function for continuous data
tf_Y_params = {'meth','FFT','param',2,'mtmin',3,'subinds',-3};% test function params

maxits_wendy = 5;
WENDy_args = {'maxits_wendy',maxits_wendy,...
    'lambdas',10.^linspace(-4,0,50),'alpha',0.01,...
    'ittol',10^-4,'diag_reg',10^-6,'verbose',0};
autowendy = 0.95; % increment library approximate confidence interval with this confidence level 
X_var = [];
tol = 5; % default heuristic increment, chosen when autowendy = 0.5;
tol_min = 0.1; % lower bound on rel. residual to increment library, in case covariance severely underestimated
tol_dd_learn = 10^-8;% ODE tolerance for forward solves in computing Y(T)

pmax_IC = 4; % max poly degree for IC solve
polys_Y_Yeq = 0:3; % Y library for Yeq solve
pmax_X_Yeq = 4; % max poly degree for X terms in Yeq solve
polys_X_Xeq = 0:2; % X library in Xeq solve
pmax_Y_Xeq = 4; % max poly degree for Y terms in Xeq solve
neg_Y = 0; % toggle use negative powers for X terms in Yeq
neg_X = 0; % toggle use negative powers for Y terms in Xeq
boolT = @(T)all([min(T,[],2)>=-2 sum(T,2)>=-2 max(T,[],2)<4],2); % restrict poly terms in Yeq

dr = '/home/danielmessenger/Dropbox/Boulder/research/data/dukic collab/';
% dr = '/projects/dame8201/datasets/dukic_collab/';

load([dr,'Gregs_mod_V=0.5.mat'],'Ycell','X','t_epi','custom_tags_X',...
    'yearlength','custom_tags_Y','linregargs_fun_IC','linregargs_fun_Y',...
    'linregargs_fun_X','nstates_X','nstates_Y','W_IC_true','tags_IC_true',...
    'W_Y_true','tags_X_true','tags_Y_true','W_X_true','tags_Ext_X_true','tags_Ext_Y_true',...
    'rhs_IC_true','rhs_Y_true','rhs_X_true','tn','t','Y');
[~,I] = findpeaks(X(:,1));

for train_time_frac = [0.75] %<<< sweep over
    if train_time_frac == 0.5
        subsamp_ts = [1 2];
    elseif train_time_frac == 0.75
        subsamp_ts = [2];
        snr_Ys = [0.01 0.05];
    elseif train_time_frac == 1
        subsamp_ts = [1 2 4 6];
        snr_Ys = [0 0.005 0.01 0.05];
    end
for subsamp_t = subsamp_ts %<<< sweep over
for kk=1:length(snr_Ys)
    snr_Y = snr_Ys(kk);

    results_cell = cell(length(ntrain_inds),length(rngs));
    sim_cell = cell(length(ntrain_inds),length(rngs));
    maps_cell = cell(length(ntrain_inds),length(rngs));
    coeffs_cell = cell(length(ntrain_inds),length(rngs));
    libs_cell = cell(length(ntrain_inds),length(rngs));
    
    for ii=1:length(ntrain_inds)
        parfor jj=1:length(rngs)
            disp([subsamp_t kk ii jj])
            num_train_inds = ntrain_inds(ii);      

            if num_train_inds<0
                gensamp_seed = unique(cell2mat(arrayfun(@(i) [i-2:i+1],I(1:-num_train_inds)','uni',0)));
            else
                gensamp_seed = rngs(jj);
            end
            noise_seed = rngs(jj);

            %%% get data
            [Y_train,X_train,train_inds,train_time,nstates_X,nstates_Y,X_in,sigma_X,sigma_Y,nX,nY] = ...
                format_data(Ycell,X,t_epi,subsamp_t,train_time_frac,num_train_inds,...
                test_length,snr_X,snr_Y,noise_alg_X,noise_alg_Y,gensamp_seed,noise_seed);
    
            tic,
            %%% run alg
            [rhs_IC,W_IC,rhs_Y,W_Y,rhs_X,W_X,lib_Y_IC,lib_X_IC,lib_Y_Yeq,lib_X_Yeq,lib_Y_Xeq,lib_X_Xeq] = ...
                wsindy_eco_fcn(toggle_zero_crossing,stop_tol,phifun_Y,tf_Y_params,WENDy_args,autowendy,tol,tol_min,tol_dd_learn,pmax_IC,polys_Y_Yeq,polys_X_Xeq,pmax_X_Yeq,pmax_Y_Xeq,neg_Y,neg_X,boolT,...
                Y_train,X_train,X_var,train_inds,train_time,t_epi,yearlength,nstates_X,nstates_Y,X_in,nX,nY,...
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

            n_err_tol = zeros(1,num_sim+1);
            x0s = [X(1,:);mean(X_train.*nX).*(1 + sqrt(3)*oos_std*(rand(num_sim,2)-0.5)*2)];
            X_pred_temp = cell(size(x0s,1),1);
            X_test_temp = cell(size(x0s,1),1);
            for j=1:size(x0s,1)
                x0 = x0s(j,:);
                num_gen = floor(size(X,1));
                sig_tmax = 0;
                num_t_epi = 112;
                if isequal(x0,X(1,:))
                    X_test = X; Ycell_test = Ycell; Y_test = Y; t_epi_test = t_epi; tn_test = tn; t_test = t;
                else
                    [X_test,Ycell_test,Y_test,t_test,tn_test,t_epi_test] = sim_hybrid_fcn(rhs_IC_true,rhs_Y_true,rhs_X_true,x0,nstates_Y,...
                    num_gen,num_t_epi,yearlength,sig_tmax,tol_dd_sim,toggle_zero_crossing,inf);
                end
                [X_pred,Ycell_pred,Y_pred,t_pred,tn_pred,t_epi_pred] = sim_hybrid_fcn(rhs_IC,rhs_Y,rhs_X,x0,nstates_Y,...
                    num_gen,num_t_epi,yearlength,sig_tmax,tol_dd_sim,toggle_zero_crossing,stop_tol*max(max(cell2mat(Y_train)./nY)));
            
                n = size(X_pred,1);
                try
                    cumerr = arrayfun(@(i)norm(vecnorm(X_pred(1:i,:)-X_test(1:i,:),2,2))/norm(vecnorm(X_test(1:i,:),2,2)),(1:n)');
                catch
                    cumerr = err_tol+1;
                end
                n_err_tol_temp = find(cumerr>err_tol,1);
                if isempty(n_err_tol_temp)
                    n_err_tol(j) = n-1;
                else
                    n_err_tol(j) = n_err_tol_temp-1;
                end
                X_pred_temp{j} = X_pred;
                X_test_temp{j} = X_test;
            end

            results_cell{ii,jj} = [errs_2_IC,errs_inf_IC,tpr_IC,...
                errs_2_Y,errs_inf_Y,tpr_Y,...
                errs_2_X,errs_inf_X,tpr_X,RT,n_err_tol];

            sim_cell{ii,jj} = {X_test_temp,X_pred_temp};
            maps_cell{ii,jj} = {rhs_IC,rhs_Y,rhs_X};
            coeffs_cell{ii,jj} = {W_IC,W_Y,W_X};
            libs_cell{ii,jj} = cellfun(@(L)cell2mat(L.tags'),{lib_Y_IC,lib_X_IC,lib_Y_Yeq,lib_X_Yeq,lib_Y_Xeq,lib_X_Xeq},'uni',0);
    
        end
    end
    if all(ntrain_inds<0)
        save([dr,'sweep_snrY_',num2str(snr_Y),'_ttf_',num2str(train_time_frac),'_subt_',num2str(subsamp_t),'_mits_',num2str(maxits_wendy),'_peaks.mat'])
    else
        save([dr,'sweep_snrY_',num2str(snr_Y),'_ttf_',num2str(train_time_frac),'_subt_',num2str(subsamp_t),'_mits_',num2str(maxits_wendy),'.mat'])
    end
end
end
end

function [value, isterminal, direction] = myEvent(T, Y, thresh, toggle_zero_crossing)
    if toggle_zero_crossing
        value      = or(norm(Y) >= thresh, any(Y==0));
    else
        value      = norm(Y) >= thresh;
    end
    isterminal = 1;   % Stop the integration
    direction  = 0;
end