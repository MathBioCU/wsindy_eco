%% feedback FitzHugh-Nagumo

nstates_Y = 2;
nstates_X = 1;
nY = ones(1,nstates_Y);
nX = ones(1,nstates_X);
tags_IC_true = [0;1;2];
W_IC_true = {[0.1 -0.5 0],[0 0 1]};
R = -0.1;
tags_X_true = [0;1];
tags_Y_true = [[1 0];[3 0];[0 1];[0 0]] ;
W_Y_true = {[[3 0];[-3 0];[3 0];[0 R]],[[-1/3 0];[0 0];[1/15 0];[17/150 0]]};

tags_Ext_X_true = [0;1;2];
tags_Ext_Y_true = [[0 0];[1 0];[0 1]];
W_X_true = {[[1 0 0];[0 1 -3];[0 0 -3]]};

tags_Ext_X_true = [1;2];
tags_Ext_Y_true = [[0 0];[1 0]];
W_X_true = {[[3.7 0.2];[-3.7 0.2]]};

[rhs_IC_true,W_IC,tags_X_IC] = shorttime_map(W_IC_true,zeros(1,nstates_Y),tags_IC_true,nY,nX);
rhs_IC_true = @(X) rhs_IC_true(zeros(nstates_Y,1),X(:));
[rhs_Y_true,W_Y,tags_Y_Yeq] = shorttime_map(W_Y_true,tags_Y_true,tags_X_true,nY,nX);
[rhs_X_true,W_X] = longtime_map(W_X_true,tags_Ext_X_true,tags_Ext_Y_true,nX,nY);

tmax_epi = 8;
x0 = 0.1;
num_gen = 30;
Ycell = cell(num_gen-1,1);
t_epi = cell(num_gen-1,1);
X = zeros(num_gen,nstates_X);
yearlength = tmax_epi;
sig_tmax = 0;
X(1,:) = x0;
num_t_epi = 256;
tol_dd = 10^-10;
options_ode_sim = odeset('RelTol',tol_dd,'AbsTol',tol_dd*ones(1,nstates_Y));
for n=1:num_gen-1
    tmax_epi_rand = tmax_epi + 2*(rand-0.5)*sig_tmax;
    t_epi{n} = linspace(0,tmax_epi_rand,num_t_epi)';
    rhs_Y = @(y)rhs_Y_true(y,X(n,:));
    x0 = rhs_IC_true(X(n,:));
    [~,x] = ode15s(@(t,x)rhs_Y(x),t_epi{n},x0,options_ode_sim);
    X(n+1,:) = rhs_X_true(x(end,:),X(n,:));
    Ycell{n} = x;
end
tn = (0:num_gen-1)*yearlength;
Y = cell2mat(Ycell); 
t = cell2mat(arrayfun(@(i)(i-1)*yearlength+t_epi{i},(1:num_gen-1)','uni',0));
for j=1:nstates_Y
    subplot(nstates_Y,1,j)
        plot(tn, X,'b-o',t, Y(:,j),'r-','linewidth',3)
end

custom_tags_Y = []; %<<< decision needs to made
custom_tags_X = []; %<<< decision needs to made
linregargs_fun_IC = @(G,b){};%{'Aineq',-G,'bineq',zeros(size(G,1),1)};
linregargs_fun_Y = @(G,b){};
linregargs_fun_X = @(G,b){};%{'Aineq',-G,'bineq',zeros(size(G,1),1)};

function [rhs_X,w_X] = longtime_map(w_X,tags_Ext_NZ,tags_Ext_Y,nX,nY)
    supp_X = cellfun(@(W)find(any(W~=0,2)),w_X(:),'uni',0);
    features_Y = cell(length(w_X),1);
    for i=1:length(w_X)
        for j=1:length(find(supp_X{i}))
            f = @(Y) Y(1)*0;
            for k=1:size(tags_Ext_Y,1)
                  w_X{i}(supp_X{i}(j),k) = [w_X{i}(supp_X{i}(j),k)/prod(nX.^tags_Ext_NZ(supp_X{i}(j),:))*nX(i)]*prod((1./nY).^tags_Ext_Y(k,:));
                  f = @(Y) f(Y)+w_X{i}(supp_X{i}(j),k)*prod(Y.^tags_Ext_Y(k,:));
            end
            features_Y{i}{j} = f;
        end
    end    
    features_NZ = cellfun(@(s)arrayfun(@(i)@(X)prod(X.^tags_Ext_NZ(i,:)),s,'uni',0),supp_X,'uni',0);
    param_map_X = @(Y)cellfun(@(f)arrayfun(@(i)f{i}(Y),1:length(f)),features_Y,'uni',0);
    rhs_X = @(Y,X) max(rhs_fun_1inp(features_NZ,param_map_X(Y),X),0);
end
function dx = rhs_fun_1inp(features,params,x)
    dx = zeros(length(params),1);
    for i=1:length(params)
        if ~isempty(features{i})
            dx(i) = dot(cellfun(@(z1) z1(x),features{i}),params{i});
        end
    end
end