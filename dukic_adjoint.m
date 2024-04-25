%%% note: throughout assumes same libraries for each component

%%% sim params
n = 1;
x = X_train(X_in(n),:);
y_data = Y_train{n};
tgrid = train_time;
tol_dd = 10^-6;
S_IC = WS_IC.get_supp;
num_IC = length(cell2mat(S_IC));
S_Y = WS_Yeq.get_supp;
% w0 = [WS_IC.weights(WS_IC.weights~=0);WS_Yeq.weights(WS_Yeq.weights~=0)];

w = w0;
grad_step = 0.1;
numits = 1000;
for j=1:numits
    %%% get y_sim
    WS_IC.weights(WS_IC.weights~=0) = w(1:num_IC);
    W_IC = cellfun(@(w)reshape(w,length(lib_X_IC.terms),[])',WS_IC.reshape_w,'uni',0);
    rhs_IC = shorttime_map(W_IC,lib_Y_IC,lib_X_IC,nX,nY);
    rhs_IC = @(X) rhs_IC(zeros(nstates_Y,1),X(:));
    WS_Yeq.weights(WS_Yeq.weights~=0) = w(num_IC+1:end);
    W_Y = cellfun(@(w)reshape(w,length(lib_X_Yeq.terms),[])',WS_Yeq.reshape_w,'uni',0);
    rhs_Y = shorttime_map(W_Y,lib_Y_Yeq,lib_X_Yeq,nX,nY);
    options_ode_sim = odeset('RelTol',tol_dd_learn,'AbsTol',tol_dd_learn*ones(1,nstates_Y),'Events',@(T,Y)myEvent(T,Y,5*max(max(cell2mat(Y_train)))));
    rhs_learned = @(y)rhs_Y(y,x.*nX);
    Y0 = rhs_IC(x.*nX);
    [t_n,y_sim,TE]=ode15s(@(t,x)rhs_learned(x),tgrid,Y0,options_ode_sim);

    [funval(j),gradY] = fun(w,tgrid,y_sim./nY,y_data,x,WS_Yeq,WS_IC,tol_dd);
    for k=1:nstates_Y
        subplot(nstates_Y,1,k)
        plot(tgrid,y_sim(:,k)./nY(k),tgrid,y_data(:,k))
        drawnow
    end
    funval(j)
    % [w grad_step*gradY]
    w = w + grad_step*gradY;
end



function [funval,gradfun] = fun(w,tgrid,y_sim,y_data,x,WS_Yeq,WS_IC,tol_dd)
    funval = mean(vecnorm(y_sim-y_data,2,2));

    %%% gather vectors
    S_IC = WS_IC.get_supp;
    num_IC = length(cell2mat(S_IC));
    S_Y = WS_Yeq.get_supp;
    
    WS_Yeq.weights(WS_Yeq.weights~=0) = w(num_IC+1:end);
    w_Y_sparse = cellfun(@(w,s)w(s),WS_Yeq.reshape_w,S_Y,'uni',0);

    %%% initialize adjoint
    nstates_Y = size(y_sim,2);
    P0 = zeros(nstates_Y,num_IC+length(find(WS_Yeq.weights)));
    ind = 0;
    for j=1:nstates_Y
        P0(j,ind+(1:length(S_IC{j}))) = WS_IC.lib(j).evalterms([y_sim(1,:) x],S_IC{j});
        ind = ind + length(S_IC{j});
    end

    RHS = @(t,P) rhs(t,tgrid,y_sim,x,P,WS_Yeq,num_IC,S_Y,w_Y_sparse);
    options_ode_sim = odeset('RelTol',tol_dd,'AbsTol',tol_dd*ones(1,length(P0(:))));

    [~,P_n]=ode15s(RHS,tgrid,reshape(P0',[],1),options_ode_sim);
    gradfun = arrayfun(@(i) reshape(P_n(i,:)',[],nstates_Y)*(y_sim(i,:)-y_data(i,:))',1:length(tgrid),'uni',0);
    gradfun = mean(cell2mat(gradfun),2);
end

function Pdot = rhs(t,tgrid,y,x,P,WS_Yeq,num_IC,S_Y,w_Y_sparse)
    y_val = interp1(tgrid,y,t);
    delYdot = get_delYdot(y_val,x,WS_Yeq,num_IC,S_Y,w_Y_sparse);
    linprop = get_linprop(y_val,x,WS_Yeq,S_Y,w_Y_sparse);
    Pdot = linprop*reshape(P,length(y_val),[]);
    Pdot = Pdot(:)+delYdot;
end

%%% get forcing
function delYdot = get_delYdot(y,x,WS_Yeq,num_IC,S_Y,w_Y_sparse)
    nstates_Y = length(y);
    delYdot = zeros(nstates_Y,num_IC+length(find(WS_Yeq.weights)));
    ind = num_IC;
    for j=1:nstates_Y
        % [y(:)' x(:)']
        A = WS_Yeq.lib(j).evalterms([y(:)' x(:)'],S_Y{j});
        delYdot(j,ind+(1:length(w_Y_sparse{j}))) = A.*w_Y_sparse{j}';
        ind = ind + length(w_Y_sparse{j});
    end
    delYdot = reshape(delYdot',[],1);
end

%%% get linear propagator
function linprop = get_linprop(y,x,WS_Yeq,S_Y,w_Y_sparse)
    nstates_Y = length(y);
    linprop = zeros(nstates_Y);
    for j=1:nstates_Y
        A = WS_Yeq.lib(j).evalGradterms([y(:)' x(:)'],S_Y{j});
        A = cell2mat(cellfun(@(g)[g{1:nstates_Y}]',A,'uni',0)); 
        linprop(j,:) = A*w_Y_sparse{j};
    end
end