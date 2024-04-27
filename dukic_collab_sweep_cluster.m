addpath(genpath('wsindy_obj_base'))
ntrain_inds = 9:3:39;
rngs = 1:500;

snr_X = 0; %<<< sweep over 

noise_alg_X = 'logn'; noise_alg_Y = 'logn'; %<<< fixed
test_length = 40; %<<< fixed
err_tol = 0.2; %<<< fixed
stop_tol = 10; %<<< fixed
toggle_zero_crossing = 1; %<<< fixed

tol_dd_sim = 10^-10; %<<< doesn't affect alg
% phifun_Y = optTFcos(3,0);%@(t)exp(-9./(1-t.^2));%<<< user choice
phifun_Y = @(t)(1-t.^2).^9;
num_sim = 5;

% tf_Y_params = {'meth','timefrac','param',0.15,'mtmin',5,'subinds',-3};%<<< user choice
tf_Y_params = {'meth','FFT','param',2,'mtmin',3,'subinds',-3};%<<< user choice

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

% dr = '/home/danielmessenger/Dropbox/Boulder/research/data/dukic collab/';
dr = '/projects/dame8201/datasets/dukic_collab/';

load([dr,'Gregs_mod_V=0.5.mat'],'Ycell','X','t_epi','custom_tags_X',...
    'yearlength','custom_tags_Y','linregargs_fun_IC','linregargs_fun_Y',...
    'linregargs_fun_X','nstates_X','nstates_Y','W_IC_true','tags_IC_true',...
    'W_Y_true','tags_X_true','tags_Y_true','W_X_true','tags_Ext_X_true','tags_Ext_Y_true',...
    'rhs_IC_true','rhs_Y_true','rhs_X_true','tn','t','Y');

results_cell = cell(length(ntrain_inds),length(rngs));
for train_time_frac = [0.75] %<<< sweep over
    if train_time_frac == 0.5
        subsamp_ts = [1 2];
    elseif train_time_frac == 0.75
        subsamp_ts = [1];
        snr_Ys = [0.005 0.02 0.03];
    elseif train_time_frac == 1
        subsamp_ts = [1 2 4 6];
        snr_Ys = [0 0.005 0.01 0.05];
    end
for subsamp_t = subsamp_ts %<<< sweep over
for kk=1:length(snr_Ys)
    snr_Y = snr_Ys(kk);
    for ii=1:length(ntrain_inds)
        parfor jj=1:length(rngs)
            % if snr_Y==0
                gensamp_seed = rngs(jj);
            % else
                % gensamp_seed = 2024;
            % end
             disp([subsamp_t kk ii jj])
             rng(rngs(jj)); rng_seed = rng().Seed; rng(rng_seed);
    
            num_train_inds = ntrain_inds(ii);      
            %%% get data
            [Y_train,X_train,train_inds,train_time,nstates_X,nstates_Y,X_in,sigma_X,sigma_Y,nX,nY] = ...
                format_data(Ycell,X,t_epi,subsamp_t,train_time_frac,num_train_inds,...
                test_length,snr_X,snr_Y,noise_alg_X,noise_alg_Y,gensamp_seed,rng_seed);
    
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

            n_err_tol = zeros(1,num_sim+1);
            x0s = [X(1,:);mean(X_train.*nX).*(1 + sqrt(3)*0.2*(rand(num_sim,2)-0.5)*2)];
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
        
            end

            results_cell{ii,jj} = [errs_2_IC,errs_inf_IC,tpr_IC,...
                errs_2_Y,errs_inf_Y,tpr_Y,...
                errs_2_X,errs_inf_X,tpr_X,RT,n_err_tol];
            % sim_cell{ii,jj} = {X_pred,Ycell_pred};
    
        end
    end
    save([dr,'sweep_snrY_',num2str(snr_Y),'_ttf_',num2str(train_time_frac),'_subt_',num2str(subsamp_t),'.mat'])
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