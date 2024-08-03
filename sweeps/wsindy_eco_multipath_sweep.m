% ntrain_inds = [-5];
% peak_width = 2;
% rngs = 1:100;
% snr_Ys = [0.002 0.01 0.02 0.05];
% snr_X = 0;
% train_time_frac = 0.75;
% subsamp_t = 3;
% maxits_wendy = 5;

addpath(genpath('../wsindy_obj_base'))
addpath(genpath('../utils'))

noise_alg_X = 'logn'; % noise distribution for X
noise_alg_Y = 'logn'; % noise distribution for Y

test_length = 40; % number of generations to test over
err_tol = 0.5; % tol for n_tol= number of generations for which cumulative rel err < tol
stop_tol = 100; % halt simulations if values exceed max observed by this multitude
toggle_zero_crossing = 0; % halt simulations that are non-positive
toggle_scale = 1;

tol_dd_sim = 10^-12; % ODE tolerance (abs,rel) for diagnostic sim
num_sim = 5; % number of out-of-sample testing simulations
oos_std = 0.5; % std of out-of-sample ICs, uniformly randomly sampled around training IC

phifun_Y = @(t)(1-t.^2).^9; % test function for continuous data
tf_Y_params = {'meth','FFT','param',1,'mtmin',3,'subinds',-3};% test function params

WENDy_args = {'maxits_wendy',maxits_wendy,...
    'lambdas',10.^linspace(-3,-1,40),'alpha',0.01,...
    'ittol',10^-4,'diag_reg',10^-6,'verbose',0};
autowendy = 0.95; % increment library approximate confidence interval with this confidence level 
toggle_X_var = 'true';
tol = 5; % default heuristic increment, chosen when autowendy = 0.5;
tol_min = 0.1; % lower bound on rel. residual to increment library, in case covariance severely underestimated
tol_dd_learn = 10^-10;% ODE tolerance for forward solves in computing Y(T)

pmax_IC = 3; % max poly degree for IC solve
polys_Y_Yeq = 0:3; % Y library for Yeq solve
pmax_X_Yeq = 3; % max poly degree for X terms in Yeq solve
polys_X_Xeq = 0:2; % X library in Xeq solve
pmax_Y_Xeq = 3; % max poly degree for Y terms in Xeq solve
neg_Y = 0; % toggle use negative powers for X terms in Yeq
neg_X = 0; % toggle use negative powers for Y terms in Xeq
boolT = @(T)all([max(T,[],2)<=2 max(T(:,6:end),[],2)==0],2); % restrict poly terms in Yeq
boolTL = repmat({@(T,L) min(T(:,min(find(L.ftag),end)),[],2)>0},1,6);
custom_tags_Y = {[],[],[],[],[],[1 1 0 1 0 0]}; % custom Y tags for Yeq. Example: {[1+V zeros(1,nstates_Y-2) 1]}. Stored with data
custom_tags_X = {}; % custom X tags for Yeq. Stored with data
linregargs_fun_IC = @(WS){}; % addition linear regression args, including constraints, as function of WSINDy model object
linregargs_fun_Y = @(WS){}; % Stored with data
linregargs_fun_X = @(WS){}; % Stored with data

% dr = '~/Desktop/';
dr = '/projects/dame8201/datasets/dukic_collab/';

load([dr,'host_multipath_6-3_c.mat'],'Ycell','X','t_epi','yearlength',...
    'W_IC_true','W_X_true','W_Y_true','tags_Ext_X_true','tags_Ext_Y_true','tags_IC_true','tags_X_true','tags_Y_true',...
    'rhs_Y_true','rhs_X_true','rhs_IC_true') %%% load in data include true model to benchmark

num_gen = 50;

num_gen = min(num_gen,size(X,1));
X = X(1:num_gen,:);
Ycell = Ycell(1:num_gen-1);
t_epi = t_epi(1:num_gen-1);
tn = (0:num_gen-1)*yearlength; % discrete time
num_t_epi = length(t_epi{1});

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
                gensamp_seed = peak_width;
            else
                gensamp_seed = rngs(jj);
            end
            noise_seed = rngs(jj);

            %%% get data
            [Y_train,X_train,train_inds,train_time,nstates_X,nstates_Y,X_in,sigma_X,sigma_Y,nX,nY] = ...
                format_data(Ycell,X,t_epi,subsamp_t,train_time_frac,num_train_inds,...
                test_length,snr_X,snr_Y,noise_alg_X,noise_alg_Y,gensamp_seed,noise_seed,toggle_scale);
            if isequal(toggle_X_var,'true')
                X_var = max(sigma_X,0);
            else
                X_var = [];
            end
            tic,
            
            %%% run alg
            [rhs_IC,W_IC,rhs_Y,W_Y,rhs_X,W_X,lib_Y_IC,lib_X_IC,lib_Y_Yeq,lib_X_Yeq,lib_Y_Xeq,lib_X_Xeq] = ...
                wsindy_eco_fcn(toggle_zero_crossing,stop_tol,phifun_Y,tf_Y_params,WENDy_args,autowendy,tol,tol_min,...
                tol_dd_learn,pmax_IC,polys_Y_Yeq,polys_X_Xeq,pmax_X_Yeq,pmax_Y_Xeq,neg_Y,neg_X,boolT,boolTL,...
                Y_train,X_train,X_var,train_inds,train_time,t_epi,yearlength,nstates_X,nstates_Y,X_in,nX,nY,...
                custom_tags_X,custom_tags_Y,linregargs_fun_IC,linregargs_fun_Y,linregargs_fun_X);
            RT = toc;

            %%% process results
            W_IC_compare = inject_coeff_param(W_IC_true(:),zeros(1,nstates_Y),tags_IC_true,cell2mat(lib_Y_IC.tags'),cell2mat(lib_X_IC.tags'));
            errs_2_IC = cellfun(@(W,V) norm(W-V)/norm(V),W_IC,W_IC_compare);
            errs_inf_IC = cellfun(@(W,V) max(abs(W(V~=0)-V(V~=0))./abs(V(V~=0))),W_IC,W_IC_compare,'un',0);
            errs_inf_IC = cell2mat(errs_inf_IC(cellfun(@(e)~isempty(e),errs_inf_IC)));
            tpr_IC = cellfun(@(W,V)tpscore(W,V),W_IC,W_IC_compare);
            
            W_Y_compare = arrayfun(@(i)...
                cell2mat(inject_coeff_param(W_Y_true(i),tags_Y_true,tags_X_true, cell2mat(lib_Y_Yeq(i).tags'),cell2mat(lib_X_Yeq.tags'))),...
                (1:length(W_Y_true))','un',0);
            errs_2_Y = cellfun(@(W,V) norm(W-V)/norm(V),W_Y,W_Y_compare);
            errs_inf_Y = cellfun(@(W,V) max(abs(W(V~=0)-V(V~=0))./abs(V(V~=0))),W_Y,W_Y_compare,'un',0);
            errs_inf_Y = cell2mat(errs_inf_Y(cellfun(@(e)~isempty(e),errs_inf_Y)));
            tpr_Y = cellfun(@(W,V)tpscore(W,V),W_Y,W_Y_compare);
            
            W_X_compare = inject_coeff_param(W_X_true(:),tags_Ext_X_true,tags_Ext_Y_true,cell2mat(lib_X_Xeq.tags'),cell2mat(lib_Y_Xeq.tags'));
            errs_2_X = cellfun(@(W,V) norm(W-V)/norm(V),W_X,W_X_compare);
            errs_inf_X = cellfun(@(W,V) max(abs(W(V~=0)-V(V~=0))./abs(V(V~=0))),W_X,W_X_compare,'un',0);
            errs_inf_X = cell2mat(errs_inf_X(cellfun(@(e)~isempty(e),errs_inf_X)));
            tpr_X = cellfun(@(W,V)tpscore(W,V),W_X,W_X_compare);

            n_err_tol = zeros(1,num_sim+1);
            % x0s = [X(1,:);mean(X_train.*nX).*(1 + sqrt(3)*oos_std*(rand(num_sim,2)-0.5)*2)];

            x0s = [X(train_inds(1),:);X(min(train_inds(1)+1,end),:);mean(X_train.*nX).*(1 + sqrt(3)*oos_std*(rand(num_sim,nstates_X)-0.5)*2)];
  
            X_pred_temp = cell(size(x0s,1),1);
            X_test_temp = cell(size(x0s,1),1);
            for j=1:size(x0s,1)
                disp([subsamp_t kk ii jj j])
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

            results_cell{ii,jj} = {errs_2_IC,errs_inf_IC,tpr_IC,...
                errs_2_Y,errs_inf_Y,tpr_Y,...
                errs_2_X,errs_inf_X,tpr_X,RT,n_err_tol};

            sim_cell{ii,jj} = {X_test_temp,X_pred_temp};
            maps_cell{ii,jj} = {rhs_IC,rhs_Y,rhs_X};
            coeffs_cell{ii,jj} = {W_IC,W_Y,W_X};
            libs_cell{ii,jj} = cellfun(@(L)cell2mat(L.tags'),{lib_Y_IC,lib_X_IC,lib_X_Yeq,lib_Y_Xeq},'uni',0);
    
        end
    end
    if all(ntrain_inds<0)
        save([dr,'TwoPath_snrY_',num2str(snr_Y),'_ttf_',num2str(train_time_frac),'_subt_',num2str(subsamp_t),'_mits_',num2str(maxits_wendy),'_peaks.mat'])
    else
        save([dr,'TwoPath_snrY_',num2str(snr_Y),'_ttf_',num2str(train_time_frac),'_subt_',num2str(subsamp_t),'_mits_',num2str(maxits_wendy),'.mat'])
    end
end