%% get samples
n= 200;
X_pred_cloud = cell(n,1);
test_ind = 1;
x0 = X_test_cell{test_ind}(1,:);
% x0s = [X(train_inds(1),:);mean(X_train.*nX).*(1 + sqrt(3)*oos_std*(rand(num_sim,2)-0.5)*2)];

parfor j=1:n
    rng('shuffle')
    disp(j)
    W_IC_s = wendy_hybrid_sample(W_IC,CovW_IC);
    W_Y_s = wendy_hybrid_sample(W_Y,CovW_Y);
    W_X_s = wendy_hybrid_sample(W_X,CovW_X);
    
    rhs_IC_s = shorttime_map(W_IC_s,lib_Y_IC,lib_X_IC,ones(1,nstates_Y),ones(1,nstates_X));
    rhs_IC_s = @(X) rhs_IC_s(zeros(nstates_Y,1),X(:));
    rhs_Y_s = shorttime_map(W_Y_s,lib_Y_Yeq,lib_X_Yeq,ones(1,nstates_Y),ones(1,nstates_X));
    rhs_X_s = longtime_map(W_X_s,lib_X_Xeq,lib_Y_Xeq,ones(1,nstates_X),ones(1,nstates_Y));
    
    [X_pred_cloud{j},Ycell_pred,Y_pred,t_pred,tn_pred,t_epi_pred] = sim_hybrid_fcn(...
        rhs_IC_s,rhs_Y_s,rhs_X_s,x0,...
        nstates_Y,num_gen,num_t_epi,yearlength,sig_tmax,...
        tol_dd_sim,toggle_zero_crossing,...
        stop_tol*max(max(cell2mat(Y_train)./nY)));    
end

%% view prediction distributions

ind = 1;
X_ext = cellfun(@(X)[X;NaN*ones(num_gen-size(X,1),nstates_X)],X_pred_cloud,'Un',0);
clf
Nmin = 10;
Nmax = num_gen;% num_gen;%size(X_pred,1);
gap = 3;
toggle_same_height = 0;

figure(1);clf
hold on
cellfun(@(x)plot(log10(x(:,ind)),0:size(x,1)-1,'linewidth',0.5,'color',[0.8 0.8 0.8]),X_pred_cloud)
for j=Nmax:-gap:Nmin
    if j~=1 % don't plot IC, will be delta
        dat = cellfun(@(X)X(j,ind),X_ext);
        dat = dat(dat>=0);
        dat = log10(dat);
        [h,xi] = ksdensity(dat(~isnan(dat)));
        if toggle_same_height==1
            hm = gap/max(h);
        else
            if j==Nmax
                hm = gap/max(h)/5;
            end
        end
        fill(xi,h*hm+j-1,[1 0.6 0.3],'linewidth',2) ;
        hold on
    end
end
hh=plot(log10(X_test_cell{test_ind}(:,ind)),0:num_gen-1,'go-',log10(X_pred_cell{test_ind}(:,ind)),0:size(X_pred_cell{test_ind},1)-1,'b--','linewidth',4,'markersize',12);
legend(hh,{'true model output','learned model output'},'location','sw','interpreter','latex');
ylim([Nmin-1 Nmax])
view([90 -90])

Xall = cell2mat(X_pred_cloud); Xall = Xall(:,ind);
xr = log10(min(Xall,[],'omitnan'))-2:2:log10(max(Xall,[],'omitnan'))+2;
set(gca,'Xtick',xr,'Xticklabels',num2str(10.^xr','%0.0e'))
if ind == 1
    xlabel('Host density ($N_n$)','interpreter','latex')
else
    xlabel('Path. density ($Z_n$)','interpreter','latex')
end
ylabel('Generation ($n$)','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','fontsize',18)
set(gcf,'position',[1970         232        1249         689])