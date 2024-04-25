%%% get data
dr = '/home/danielmessenger/Dropbox/Boulder/research/data/dukic collab/';
% load([dr,'FH_feedback.mat']);
load([dr,'Gregs_mod_V=0.5.mat'],'Ycell','X','t_epi','custom_tags_X',...
    'yearlength','custom_tags_Y','linregargs_fun_IC','linregargs_fun_Y',...
    'linregargs_fun_X','W_IC_true','tags_IC_true',...
    'W_Y_true','tags_X_true','tags_Y_true','W_X_true','tags_Ext_X_true','tags_Ext_Y_true',...
    'rhs_IC_true','rhs_Y_true','rhs_X_true','tn','t','Y');

seed1 = 2024;
seed2 = 'shuffle';
snr_X = 0.0; %<<< sweep over 
snr_Y = 0.01; %<<< sweep over
train_time_frac = 1; %<<< sweep over
subsamp_t = 111;
num_train_inds = 30; %<<< sweep over

noise_alg_X = 'logn'; 
noise_alg_Y = 'logn'; %<<< fixed
test_length = 50; %<<< fixed

WA = {{'maxits_wendy',10,...
    'lambdas',10.^linspace(-4,0,50),...
    'ittol',10^-4,'diag_reg',10^-6,'verbose',0};{'maxits_wendy',0,...
    'lambdas',10.^linspace(-4,0,50),...
    'ittol',10^-4,'diag_reg',10^-6,'verbose',0}};
autowendy = 1; %<<< decision needs to made
tol = 5; %<<< decision needs to made
tol_min = 0.1; %<<< decision needs to made

polys_X_Xeq = 0:2; %<<< decision needs to made
pmax_X_Yeq = 4; %<<< decision made 
pmax_Y_Xeq = 4; %<<< decision made
neg_X = 0; %<<< decision needs to made
boolT = @(T)all([min(T,[],2)>=-2 sum(T,2)>=-2 max(T,[],2)<4],2); %<<< decision made

runs = 50;clc
% results = repmat({zeros(runs,3)},1,length(WA));
%%
for r=1:runs
[Y_train,X_train,train_inds,train_time,nstates_X,nstates_Y,X_in,sigma_X,sigma_Y,nX,nY] = ...
    format_data(Ycell,X,t_epi,subsamp_t,train_time_frac,num_train_inds,test_length,...
    snr_X,snr_Y,noise_alg_X,noise_alg_Y,seed1,seed2);
toggle_sim = size(X,1);

tic;
Yend = cell2mat(cellfun(@(y)y(end,:),Y_train,'uni',0));
E = eye(nstates_Y+nstates_X);

X_Yend = zeros(max(train_inds),nstates_X+nstates_Y);
X_Yend(train_inds,1:nstates_X) = X_train;
subinds = train_inds(X_in);
X_Yend(subinds,nstates_X+1:end) = Yend;
Uobj_X_Yend = wsindy_data(X_Yend,(0:max(train_inds)-1)*yearlength);
Xn_cov = X_Yend*0; 
sig_X = repmat(cell2mat(sigma_X),size(subinds,1),1);
sig_Y = cell2mat(cellfun(@(y)cell2mat(y),sigma_Y,'uni',0));
Xn_cov(subinds,:) = [sig_X sig_Y];
Uobj_X_Yend.R0 = spdiags(Xn_cov(:),0,numel(Xn_cov),numel(Xn_cov));

lib_X_Xeq = library('tags',get_tags(polys_X_Xeq,[],nstates_X));
lib_Y_Xeq = library();
tf_X = testfcn(Uobj_X_Yend,'meth','direct','param',1/2,'phifuns','delta','mtmin',0,'subinds',subinds);
lhs = arrayfun(@(i)term('ftag',E(i,:),'linOp',1),(1:nstates_X)','uni',0);

fprintf('\n-- run %u --\n',r)

for ii=1:length(WA)
wendy_args=WA{ii};
[rhs_X,W_X,WS_Xeq,lib_Y_Xeq,loss_X,lambda_X,w_its,res_X,res_0_X,CovW_X] = ...
    hybrid_MI(pmax_Y_Xeq,lib_X_Xeq,lib_Y_Xeq,nstates_X,nstates_Y,Uobj_X_Yend,tf_X,lhs,wendy_args,linregargs_fun_X,autowendy,tol,tol_min,nX,nY);

W_X_compare = inject_coeff_param(W_X_true,tags_Ext_X_true,tags_Ext_Y_true,cell2mat(lib_X_Xeq.tags'),cell2mat(lib_Y_Xeq.tags'));
errs_2_X = norm(reshape([W_X{:}]-[W_X_compare{:}],[],1))/norm(reshape([W_X_compare{:}],[],1));
fprintf('Xn coeff err2: %1.3e \n',errs_2_X)
errs_inf_X = abs(reshape([W_X{:}]-[W_X_compare{:}],[],1))./abs(reshape([W_X_compare{:}],[],1));
errs_inf_X = max(errs_inf_X(errs_inf_X<inf),[],'all','omitnan');
if isempty(errs_inf_X)
    errs_inf_X = NaN;
end
fprintf('Xn coeff errInf: %1.3e \n',errs_inf_X)
tpr_X = tpscore(reshape([W_X{:}],[],1),reshape([W_X_compare{:}],[],1));
fprintf('Xn TPR: %0.3f \n',tpr_X)
fprintf('\n runtime: %2.3f \n',toc)
results{ii}(r,:) = [errs_2_X errs_inf_X tpr_X];
end
end

%%
load([dr,'X_map_test.mat'],'results','snr_Y','num_train_inds')
boxplot([results{1}(:,3) results{2}(:,3)])
set(gca,'Xticklabels',{'WENDy','OLS'})
title(['TPR: ',num2str(snr_Y*100),'% noise in Y; ',num2str(num_train_inds),' observed gens'])
