% dr = '/home/danielmessenger/Dropbox/Boulder/research/data/dukic collab/';
% % load([dr,'FH_feedback.mat']);
% load([dr,'Gregs_mod_V=0.5.mat'])
% tol_dd=tol_ode; toggle_zero_crossing=0;stop_tol=inf; num_gen = M; x0 = X0s(1:2);

% [X,Ycell,Y,t,tn] = sim_hybrid_fcn(rhs_IC_true,rhs_Y_true,rhs_X_true,x0,nstates_Y,...
%     num_gen,num_t_epi,yearlength,sig_tmax,tol_dd,toggle_zero_crossing,stop_tol);
% plot(X)

function [X,Ycell,Y,t,tn,t_epi] = sim_hybrid_fcn(rhs_IC,rhs_xy,rhs_X,x0,nstates_Y,...
    num_gen,num_t_epi,yearlength,sig_tmax,tol_dd,toggle_zero_crossing,stop_tol)
    Ycell = {};
    t_epi = {};
    X = x0(:)';
    options_ode_sim = odeset('RelTol',tol_dd,'AbsTol',tol_dd*ones(1,nstates_Y),...
        'Events',@(T,Y)myEvent(T,Y,stop_tol,toggle_zero_crossing));
    n=1;
    check=1;
    t = [];
    while and(n<num_gen,check)
        tmax_epi_rand = yearlength + 2*(rand-0.5)*sig_tmax;
        t_epi = [t_epi;{linspace(0,tmax_epi_rand,num_t_epi)'}];
        rhs_Y = @(y)rhs_xy(y,X(n,:));
        [t_n,Y_n,TE] = ode15s(@(t,x)rhs_Y(x),t_epi{n},rhs_IC(X(n,:)),options_ode_sim);
        X(n+1,:) = rhs_X(X(n,:),Y_n(end,:));
        Ycell = [Ycell;{Y_n}];
        t = [t;t_n+(n-1)*yearlength];
        if toggle_zero_crossing==1
            check = all([all([isempty(TE) X(n+1,:)>=0]) ~myEvent(0,X(n+1,:),stop_tol,toggle_zero_crossing)]);
        else
            check = all([isempty(TE) ~myEvent(0,X(n+1,:),stop_tol,toggle_zero_crossing)]);
        end
        n=n+1;
    end
    tn = (0:n-1)*yearlength;
    Y = cell2mat(Ycell); 
    % t = cell2mat(arrayfun(@(i)(i-1)*yearlength+t_epi{i},(1:n-1)','uni',0));
end

function [value, isterminal, direction] = myEvent(T, Y, thresh, toggle_zero_crossing)
    if toggle_zero_crossing
        value      = or(norm(Y) >= thresh, any(Y<0));
    else
        value      = norm(Y) >= thresh;
    end
    isterminal = 1;   % Stop the integration
    direction  = 0;
end