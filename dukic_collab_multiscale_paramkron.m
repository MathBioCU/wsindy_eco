%% load hyperparams
hybrid_ID_inputs; tic;

%% run alg
[rhs_IC,W_IC,rhs_Y,W_Y,rhs_X,W_X,lib_Y_IC,lib_X_IC,...
    lib_Y_Yeq,lib_X_Yeq,lib_Y_Xeq,lib_X_Xeq,WS_IC,WS_Yeq,WS_Xeq,...
    CovW_IC,CovW_Y,CovW_X,errs_Yend]= ...
    hybrid_MI_full(toggle_zero_crossing,stop_tol,phifun_Y,tf_Y_params,WENDy_args,autowendy,tol,tol_min,tol_dd_learn,pmax_IC,polys_Y_Yeq,polys_X_Xeq,pmax_X_Yeq,pmax_Y_Xeq,neg_Y,neg_X,boolT,...
    Y_train,X_train,train_inds,train_time,t_epi,yearlength,nstates_X,nstates_Y,X_in,nX,nY,...
    custom_tags_X,custom_tags_Y,linregargs_fun_IC,linregargs_fun_Y,linregargs_fun_X);

%% coeff compare
addpath(genpath('wsindy_obj_base'))
W_IC_compare = inject_coeff_param(W_IC_true,zeros(1,nstates_Y),tags_IC_true,cell2mat(lib_Y_IC.tags'),cell2mat(lib_X_IC.tags'));
errs_2_IC = norm(reshape([W_IC{:}]-[W_IC_compare{:}],[],1))/norm(reshape([W_IC_compare{:}],[],1));
fprintf('IC coeff err 2: %1.3e \n',errs_2_IC);
errs_inf_IC = abs(reshape([W_IC{:}]-[W_IC_compare{:}],[],1))./abs(reshape([W_IC_compare{:}],[],1));
errs_inf_IC = max(errs_inf_IC(errs_inf_IC<inf),[],'all','omitnan');
if isempty(errs_inf_IC)
    errs_inf_IC = NaN;
end
fprintf('IC coeff errInf: %1.3e \n',errs_inf_IC)
tpr_IC = tpscore(reshape([W_IC{:}],[],1),reshape([W_IC_compare{:}],[],1));
fprintf('IC TPR: %0.3f \n',tpr_IC)

W_Y_compare = inject_coeff_param(W_Y_true,tags_Y_true,tags_X_true, cell2mat(lib_Y_Yeq.tags'),cell2mat(lib_X_Yeq.tags'));
errs_2_Y = norm(reshape([W_Y{:}]-[W_Y_compare{:}],[],1))/norm(reshape([W_Y_compare{:}],[],1));
fprintf('Ydot coeff err: %1.3e \n',errs_2_Y)
errs_inf_Y = abs(reshape([W_Y{:}]-[W_Y_compare{:}],[],1))./abs(reshape([W_Y_compare{:}],[],1));
errs_inf_Y = max(errs_inf_Y(errs_inf_Y<inf),[],'all','omitnan');
if isempty(errs_inf_Y)
    errs_inf_Y = NaN;
end
fprintf('Ydot coeff errInf: %1.3e \n',errs_inf_Y)
tpr_Y = tpscore(reshape([W_Y{:}],[],1),reshape([W_Y_compare{:}],[],1));
fprintf('Ydot TPR: %0.3f \n',tpr_Y)

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

%% sim full system
addpath(genpath('wsindy_obj_base'))
if toggle_sim>0
    num_sim = 0;
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
        if toggle_zero_crossing==1
            if any([X_pred(n,:)]<0)
                fprintf(['negative values detected in X(%u) \n'],find([X_pred(n,:)]<0))
            end
            if any(Ycell_pred{n-1}(:)<0)
                fprintf(['negative values detected in Y, %u \n'],find(Ycell_pred{n-1}(:)<0))
            end
        end
        
        fprintf('X err: %1.3e \n',vecnorm(X_pred(1:n,:)-X_test(1:n,:))./vecnorm(X_test(1:n,:)))
        cumerr = arrayfun(@(i)norm(vecnorm(X_pred(1:i,:)-X_test(1:i,:),2,2))/norm(vecnorm(X_test(1:i,:),2,2)),(1:n)');
        n_err_tol = find(cumerr>err_tol,1);
        if isempty(n_err_tol)
            n_err_tol = n-1;
        else
            n_err_tol = n_err_tol-1;
        end
        fprintf('n_{%0.1f}=%0.2f \n',err_tol,n_err_tol)
        
        if toggle_vis
            figure(j);clf;set(gcf,'position',[1755         156        1837         716])
            dukic_collab_sol_vis;
        end
    end
end
%% functions

function [value, isterminal, direction] = myEvent(T, Y, thresh, toggle_zero_crossing)
    if toggle_zero_crossing
        value      = or(norm(Y) >= thresh, any(Y==0));
    else
        value      = norm(Y) >= thresh;
    end
    isterminal = 1;   % Stop the integration
    direction  = 0;
end
% WS = wsindy_model(Uobj,lib_Y,tf);
% A_NZ = lib_NZ.evalterms(NZ_train);
% [W,G,b] = WS_opt().MSTLS_param(WS,A_NZ);
% err = max(cellfun(@(G,W,b)norm(G*reshape(W',[],1)-b)./norm(b),G,W',b));
% %% compare small scale model
% test_frac = 1;toggle_compare = [];%toggle_compare = find(~ismember(1:length(xcell),train_inds));if isempty(toggle_compare);toggle_compare = train_inds;end
% for n=1:length(toggle_compare)
%     options_ode_sim = odeset('RelTol',tol_dd_sim,'AbsTol',tol_dd_sim*ones(1,nstates_Y));
%     rhs_learned = @(y)rhs_xy(y,NZ(toggle_compare(n),:));
%     x0 = IC_map(NZ(toggle_compare(n),:));
%     t_train = t_epi(1:floor(end*test_frac));
%     [t_learned,x_learned]=ode15s(@(t,x)rhs_learned(x),t_train,x0,options_ode_sim);
%     figure(1);clf
%     for j=1:nstates_Y
%         subplot(nstates_Y,1,j)
%         plot(t_epi,xcell{toggle_compare(n)}(:,j),'b-',t_learned,x_learned(:,j),'--','linewidth',2)
%         U_true = xcell{toggle_compare(n)}(1:min(length(t_learned),floor(end*test_frac)),j);
%         err = norm(x_learned(:,j)-U_true)/norm(U_true);
%         title(['rel err=',num2str(err)])
%         legend({'data','learned'})
%     end
% end
% p=-1;err = tol_X*2;
% while and(err>tol_X,p<pmax_X)
%     p=p+1;
%     polys_X = 0:p;
%     tags_X_Xeq = get_tags(polys_X,[],nstates_X);
%     tags_Y_Xeq = get_tags(polys_X,[],nstates_Y);
%     Theta_X_Xeq = cell2mat(arrayfun(@(j)prod(X_n.^tags_X_Xeq(j,:),2),1:size(tags_X_Xeq),'uni',0));
%     Theta_Y_Xeq = cell2mat(arrayfun(@(j)prod(Yend.^tags_Y_Xeq(j,:),2),1:size(tags_Y_Xeq),'uni',0));    
%     Theta_X = cell2mat(arrayfun(@(i)kron(Theta_X_Xeq(i,:),Theta_Y_Xeq(i,:)),(1:size(Theta_X_Xeq,1))','uni',0));
%     W_X = arrayfun(@(i)sparsifyDynamics(Theta_X,b_X(:,i),lambda,1,0,ones(size(Theta_X,2),1),inf,1,[]),1:nstates_X,'uni',0);
%     % s = 15; ncv = 25; 
%     % w_X = WS_opt().subspacePursuitCV({repmat({A_kron},1,nstates_X),mat2cell(b,size(b,1),ones(1,size(b,2)))},s,ncv);
%     err = max(vecnorm(Theta_X*cell2mat(W_X)-b_X)./vecnorm(b_X));
%     if p==pmax_X
%         disp(['tolerance not reached'])
%     end
% end
% W_X = cellfun(@(w)reshape(w,size(Theta_Y_Xeq,2),[])',W_X,'uni',0);
% IC = cell2mat(cellfun(@(x)x(1,:),Y_train,'uni',0));
% p=-1;err = tol_IC*2;
% while and(err>tol_IC,p<pmax_IC)
%     p=p+1;
%     polys_X_IC = 0:p;
%     tags_X_IC = get_tags(polys_X_IC,[],nstates_X);
%     Theta_IC = cell2mat(arrayfun(@(j)prod(X_train(X_in,:).^tags_X_IC(j,:),2),1:size(tags_X_IC),'uni',0));
%     W_IC = arrayfun(@(i)sparsifyDynamics(Theta_IC,IC(:,i),lambda,1,0,ones(size(Theta_IC,2),1),inf,1,[]),1:nstates_Y,'uni',0);
%     W_IC = cell2mat(W_IC);
%     err = max(vecnorm(Theta_IC*W_IC-IC)./vecnorm(IC));
%     if p==pmax_IC
%         disp(['tolerance not reached'])
%     end
% end
% [rhs_IC,W_IC] = IC_func(W_IC,tags_X_IC,nY,nX);
% function wnew = inject_coeff(w,tags_old,tags_new)
%     wnew = zeros(size(tags_new,1),size(w,2));
%     for i=1:size(tags_old,1)
%         try
%             wnew(ismember(tags_new,tags_old(i,:),'rows'),:) = w(i,:);
%         catch
%             disp('term not included in library')
%             continue
%         end
%     end
% end
% 
% function [IC_map,IC_NZ] = IC_func(IC_NZ,tags_IC,nY,nX)
% nstates_Y = size(IC_NZ,2);
% for j=1:nstates_Y
%     f = @(X) X(1)*0;
%     for k=1:size(tags_IC,1)
%         IC_NZ(k,j) = nY(j)*IC_NZ(k,j)*prod((1./nX).^tags_IC(k,:));
%         f = @(X) f(X)+IC_NZ(k,j)*prod(X.^tags_IC(k,:));
%     end
%     IC_map{j} = f;
% end
% IC_map = @(X)cellfun(@(f)f(X),IC_map);
% end
% X_train = X_train.*exp(sigma_X*std(log(X_train)).*randn(size(X_train)));
% %% get wsindy_data 
% sigma_est_IC = cell(1,nstates_Y+nstates_X);
% sigma_est_YX = cell(length(X_in),1);
% sigma_est_X_Yend = cell(1,nstates_Y+nstates_X);
% if ~isempty(sigma_est_Y)
%     sigma_est_IC(1:nstates_Y) = num2cell(mean(cell2mat(cellfun(@(s)[s{:}],sigma_est_Y,'uni',0))));
%     sigma_est_YX = cellfun(@(s)[s,cell(1,nstates_X)],sigma_est_Y,'uni',0);
% else
%     sigma_est_Y = cell(length(X_in),1);
% end
% if ~isempty(sigma_est_X)
%     sigma_est_IC(nstates_Y+1:end) = sigma_est_X;
%     sigma_est_X_Yend = [sigma_est_X,cell(1,nstates_Y)];
%     sigma_est_YX = cellfun(@(s)[s(1:nstates_Y),sigma_est_X],sigma_est_YX,'uni',0);
% end
% Uobj_Y = cellfun(@(x,s)wsindy_data(x,train_time,'sigmas',s),Y_train,sigma_est_Y);
% IC = zeros(max(train_inds),nstates_Y+nstates_X);
% IC(train_inds,nstates_Y+1:end) = X_train;
% IC(train_inds(X_in),1:nstates_Y) = cell2mat(arrayfun(@(U)cellfun(@(x)x(1,:),U.Uobs),Uobj_Y,'uni',0));
% Uobj_IC = wsindy_data(IC,0:max(train_inds)-1,'sigmas',sigma_est_IC);
% if isempty(sigma_est_Y{1})
%     S = arrayfun(@(U)cell2mat(U.estimate_sigma).^2,Uobj_Y,'uni',0);
%     IC_cov = IC*0; 
%     IC_cov(train_inds(X_in),1:nstates_Y) = cell2mat(S);
%     Uobj_IC.R0 = spdiags(IC_cov(:),0,numel(IC_cov),numel(IC_cov));
% end
% Uobj_tot = arrayfun(@(i)...
%     wsindy_data([[Uobj_Y(i).Uobs{:}] repmat(X_train(X_in(i),:),Uobj_Y(i).dims,1)],train_time,...
%     'sigmas',sigma_est_YX{i}),(1:length(Uobj_Y))');
% E = eye(nstates_Y+nstates_X);
% if ~all(cellfun(@(s)isempty(s),sigma_est_X_Yend))
% sigma_est_X_Yend = cell(1,nstates_Y+nstates_X);
% sigma_est_X_Yend(nstates_X+1:end) = num2cell(mean(errs,1,'omitnan'));
% % end
% Uobj_X_Yend = wsindy_data(X_Yend,0:max(train_inds)-1,'sigmas',sigma_est_X_Yend);
    % [WS_Yeq,w] = polylibinc_w(WS_Yeq,2,0.5,100,1);
    % lib_Y_Yeq = library('tags',cell2mat(lib_Y_Yeq.tags')*w); % library
    % lib_X_Yeq = library('tags',cell2mat(lib_X_Yeq.tags')*w);     
    % WS_Yeq.cat_Gb('cat','blkdiag');
    % Ycell_pred = {};
    % X_pred = X(1,:);
    % options_ode_sim = odeset('RelTol',tol_dd_sim,'AbsTol',tol_dd_sim*ones(1,nstates_Y),'Events',@(T,Y)myEvent(T,Y,stop_tol*max(max(cell2mat(Y_train))),toggle_zero_crossing));
    % tn_pred = [];
    % n=1;
    % check=1;
    % while and(n<toggle_sim,check)
    %     rhs_learned = @(y)rhs_Y(y,X_pred(n,:));
    %     Y0 = rhs_IC(X_pred(n,:));
    %     t_train = t_epi{n};
    %     [t_n,Y_n,TE]=ode15s(@(t,x)rhs_learned(x),t_train,Y0,options_ode_sim);
    %     Ycell_pred = [Ycell_pred;{Y_n}];
    %     tn_pred = [tn_pred;t_n+(n-1)*yearlength];
    %     X_pred(n+1,:) = rhs_X(X_pred(n,:),Y_n(end,:));
    %     if toggle_zero_crossing==1
    %         check = all([isempty(TE) X_pred(n+1,:)>=0]);
    %     else
    %         check = isempty(TE);
    %     end
    %     n=n+1;
    % end
    % Y_pred = cell2mat(Ycell_pred);
