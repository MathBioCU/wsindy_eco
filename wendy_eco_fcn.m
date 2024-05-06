function [rhs_IC,W_IC,rhs_Y,W_Y,rhs_X,W_X,...
    lib_Y_IC,lib_X_IC,lib_Y_Yeq,lib_X_Yeq,lib_Y_Xeq,lib_X_Xeq,...
    WS_IC,WS_Yeq,WS_Xeq,...
    CovW_IC,CovW_Y,CovW_X,...
    errs_Yend]= ...
    wendy_eco_fcn(toggle_zero_crossing,stop_tol,phifun_Y,tf_Y_params,WENDy_args,autowendy,tol_dd_learn,...
        Y_train,X_train,X_var,train_inds,train_time,t_epi,yearlength,nstates_X,nstates_Y,X_in,nX,nY,...
        tags_IC_X,tags_Y_Yeq,tags_X_Yeq,tags_X_Xeq,tags_Y_Xeq,W_IC_sparse,W_Y_sparse,W_X_sparse)
    
    %%% get wsindy_data 
    Uobj_Y = cellfun(@(x,t)wsindy_data(x,t),Y_train,train_time);
    Uobj_tot = arrayfun(@(i)...
        wsindy_data([[Uobj_Y(i).Uobs{:}] repmat(X_train(X_in(i),:),Uobj_Y(i).dims,1)],train_time{i}),(1:length(Uobj_Y))');
    foo = arrayfun(@(U)U.estimate_sigma('set',true),Uobj_tot,'uni',0);
    for j=1:length(Uobj_tot)
        Uobj_tot(j).sigmas(nstates_Y+1:end) = num2cell(X_var(X_in(j),:));
    end
    
    IC = zeros(max(train_inds),nstates_Y+nstates_X);
    IC(train_inds,nstates_Y+1:end) = X_train;
    IC(train_inds(X_in),1:nstates_Y) = cell2mat(arrayfun(@(U)cellfun(@(x)x(1,:),U.Uobs),Uobj_Y,'uni',0));
    Uobj_IC = wsindy_data(IC,0:max(train_inds)-1);
    if autowendy>0
        S = arrayfun(@(U)cell2mat(U.estimate_sigma).^2,Uobj_Y,'uni',0);
        IC_cov = IC*0; 
        IC_cov(train_inds(X_in),1:nstates_Y) = cell2mat(S);
        IC_cov(train_inds(X_in),nstates_Y+1:end) = X_var(X_in,:).^2;
        Uobj_IC.R0 = spdiags(IC_cov(:),0,numel(IC_cov),numel(IC_cov));
    end
    E = eye(nstates_Y+nstates_X);
    
    %%% get IC_map
    linregargs = {'S',{logical(cell2mat(cellfun(@(w)w(:),W_IC_sparse(:),'Un',0)))}};
    lib_Y_IC = library('tags',zeros(1,nstates_Y));
    lib_X_IC = library('tags',tags_IC_X);
    tf_IC = testfcn(Uobj_IC,'meth','direct','param',0,'phifuns','delta','mtmin',0,'subinds',train_inds(X_in));
    lhs_IC = arrayfun(@(i)term('ftag',E(i,:)),(1:nstates_Y)','uni',0);
    [rhs_IC,W_IC,WS_IC,lib_X_IC,W_its_IC,res_IC,res_0_IC,CovW_IC] = ...
        wendy_par_fcn(lib_Y_IC,lib_X_IC,Uobj_IC,tf_IC,lhs_IC,WENDy_args,linregargs,nY,nX);
    rhs_IC = @(X) rhs_IC(zeros(nstates_Y,1),X(:));
    
    %%% get parametric small scale model
    linregargs = {'S',{logical(cell2mat(cellfun(@(w)w(:),W_Y_sparse(:),'Un',0)))}};
    tf_Yeq = arrayfun(@(U)testfcn(U,tf_Y_params{:},'phifuns',phifun_Y),Uobj_Y,'uni',0);
    lib_Y_Yeq = library('tags',tags_Y_Yeq);
    lib_X_Yeq = library('tags',tags_X_Yeq);
    lhs_Yeq = arrayfun(@(i)term('ftag',E(i,:),'linOp',1),(1:nstates_Y)','uni',0);
    [rhs_Y,W_Y,WS_Yeq,lib_X_Yeq,W_its_Y,res_Y,res_0_Y,CovW_Y] = ...
        wendy_par_fcn(lib_Y_Yeq,lib_X_Yeq,Uobj_tot,tf_Yeq,lhs_Yeq,WENDy_args,linregargs,nY,nX);
    
    %%% get Y(T)
    X_sub = find(diff(train_inds)==1);
    subinds = train_inds(X_sub);
    Yend = zeros(length(X_sub),nstates_Y);
    X_n = X_train(X_sub,:);
    errs_Yend = []; Y_new = {}; Y_ns = cell(length(X_sub),2);
    for n=1:length(X_sub)
        options_ode_sim = odeset('RelTol',tol_dd_learn,'AbsTol',tol_dd_learn*ones(1,nstates_Y),'Events',@(T,Y)myEvent(T,Y,stop_tol*max(max(cell2mat(Y_train))),toggle_zero_crossing));
        rhs_learned = @(y)rhs_Y(y,X_n(n,:).*nX);
        Y0 = rhs_IC(X_n(n,:).*nX);
        t_train = t_epi{train_inds(X_sub(n))}([1 end]);
        [t_n,Y_n,TE]=ode15s(@(t,x)rhs_learned(x),t_train,Y0,options_ode_sim);
        Y_ns(n,:) = {Y_n,t_n};
        if ismember(X_sub(n),X_in)
            try
                Y_new = [Y_new;{interp1(t_n,Y_n,train_time{X_in==X_sub(n)})./nY}];
                errs_Yend = [errs_Yend;std(Y_new{end}-Y_train{X_in==X_sub(n)})];
            catch
                Y_new = [Y_new;{NaN}];
                errs_Yend = [errs_Yend;max(abs(Y_n),'omitnan')];
            end
        end
        Yend(n,:) = Y_n(end,:);
    end
    Yend = Yend./nY;
    
    %%% get large scale model
    linregargs = {'S',{logical(cell2mat(cellfun(@(w)w(:),W_X_sparse(:),'Un',0)))}};
    X_Yend = zeros(max(train_inds),nstates_X+nstates_Y);
    X_Yend(train_inds,1:nstates_X) = X_train;
    X_Yend(subinds,nstates_X+1:end) = Yend;
    Uobj_X_Yend = wsindy_data(X_Yend,(0:max(train_inds)-1)*yearlength);
    if autowendy>0
        Xn_cov = X_Yend*0; 
        errs_Yend = fillmissing(interp1(train_inds(X_in),errs_Yend,subinds,'linear'),'linear');
        Xn_cov(subinds,1:nstates_X) = X_var(X_sub,:).^2;
        Xn_cov(subinds,nstates_X+1:end) = errs_Yend.^2;
        Uobj_X_Yend.R0 = spdiags(Xn_cov(:),0,numel(Xn_cov),numel(Xn_cov));
    end
    lib_X_Xeq = library('tags',tags_X_Xeq);
    lib_Y_Xeq = library('tags',tags_Y_Xeq);
    tf_X = testfcn(Uobj_X_Yend,'meth','direct','param',1/2,'phifuns','delta','mtmin',0,'subinds',subinds);
    lhs_X = arrayfun(@(i)term('ftag',E(i,:),'linOp',1),(1:nstates_X)','uni',0);
    [rhs_X,W_X,WS_Xeq,lib_Y_Xeq,W_its_X,res_X,res_0_X,CovW_X] = ...
        wendy_par_fcn(lib_X_Xeq,lib_Y_Xeq,Uobj_X_Yend,tf_X,lhs_X,WENDy_args,linregargs,nX,nY);
end