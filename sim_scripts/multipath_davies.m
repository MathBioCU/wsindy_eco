toggle_save = 1;
%% 2-path system
nstates_Y = 7;
nstates_X = 3;

%discrete : N,nu1,nu2,Z1,Z2
%continuos: S,nu1,nu2,P1,P2,I1,I2

yearlength = 1;
rho = -0.5;
c1 = 2.06;
c2 = 0.68;
mu1 = 0.65;%65
mu2 = 0.65;
lam = 10;
phi = 4;

nu1 = 10.5;
nu2 = 0.5;
S0 = 100;
P0 = log(lam)/nu1/yearlength;

x0 = [S0 P0 10^-6];%.*(1+max(0.2*randn(1,5),-1+eps));

tags_IC_true = [eye(3);zeros(1,3)];
W_IC_true = {[1 0 0 0],[0 1 0 0],[0 0 1 0],[0 0 0 nu1],[0 0 0 nu2],[0 0 0 0],[0 0 0 0]};

tags_X_true = [zeros(1,3);[-1 zeros(1,2)]]; 
tags_Y_true = [[1 1 0 1 0 0 0];[1 0 1 0 1 0 0];... 
                [0 1 0 2 0 0 0];[0 0 1 1 1 0 0];...
                [0 0 1 0 2 0 0];[0 1 0 1 1 0 0];...
                [0 1 0 0 0 0 0];[0 0 1 0 0 0 0]];
W_Y_true = {[[-1 -1 zeros(1,6)]',zeros(8,1)],...
            [[1 zeros(1,5) -mu1 0]',zeros(8,1)],...
            [[0 1 zeros(1,5) -mu2]',zeros(8,1)],...
            [[0 0 -c1^2 -rho*c1*c2 0 0 0 0]',zeros(8,1)],...
            [[0 0 0 0 -c2^2 -rho*c1*c2 0 0]',zeros(8,1)],...
            [zeros(8,1),[1;zeros(7,1)]],...
            [zeros(8,1),[0;1;zeros(6,1)]],...
            };

tags_Ext_X_true = [0 0 0;1 0 0];
tags_Ext_Y_true = [1 zeros(1,6);...
                    zeros(1,5) 1 0;...
                    zeros(1,6) 1;...
                    0 0 0 1 zeros(1,3);...
                    0 0 0 0 1 0 0;zeros(1,7)];
W_X_true = {[lam 0 0 0 0 0;zeros(1,6)],...
            [zeros(1,6);0 phi 0 0 0 0],...
            [zeros(1,6);0 0 phi 0 0 0],...
            };

%% sim 
addpath(genpath('../utils'));
addpath(genpath('../wsindy_obj_base'));
[rhs_IC_true,~,~] = shorttime_map(W_IC_true,zeros(1,nstates_Y),tags_IC_true,[],[]);
rhs_IC_true = @(X) rhs_IC_true(zeros(nstates_Y,1),X(:));
[rhs_Y_true,~,~] = shorttime_map(W_Y_true,tags_Y_true,tags_X_true,[],[]);
[rhs_X_true,~,~] = shorttime_map(W_X_true,tags_Ext_X_true,tags_Ext_Y_true,[],[]);

num_gen = 100;
Ycell = cell(num_gen-1,1);
t_epi = cell(num_gen-1,1);
X = zeros(num_gen,nstates_X);
sig_tmax = 0;
X(1,:) = x0;
num_t_epi = 4*yearlength*56;
options_ode_sim = odeset('RelTol',10^-10,'AbsTol',10^-10*ones(1,nstates_Y));
for n=1:num_gen-1
    tmax_epi_rand = yearlength + 2*(rand-0.5)*sig_tmax;
    t_epi{n} = linspace(0,tmax_epi_rand,num_t_epi)';
    rhs_Y = @(y)rhs_Y_true(y,X(n,:));
    [~,x] = ode15s(@(t,x)rhs_Y(x),t_epi{n},rhs_IC_true(X(n,:)),options_ode_sim);
    Ycell{n} = x;
    X(n+1,:) = rhs_X_true(X(n,:),x(end,:));
end
tn = (0:num_gen-1)*yearlength;
Y = cell2mat(Ycell); 
t = cell2mat(arrayfun(@(i)(i-1)*yearlength+t_epi{i},(1:num_gen-1)','uni',0));

figure(1)
for j=1:nstates_Y
    subplot(nstates_Y,1,j)
    semilogy(t, Y(:,j),'r-','linewidth',3)
end

figure(2)
for j=1:nstates_X
    subplot(nstates_X,1,j)
    semilogy(tn, X(:,j),'b-o','linewidth',3)
end

custom_tags_Y = [];
custom_tags_X = [-1 0 0];
linregargs_fun_IC = @(G,b){};%{'Aineq',-G,'bineq',zeros(size(G,1),1)};
linregargs_fun_Y = @(G,b){};
linregargs_fun_X = @(G,b){};

if toggle_save==1
    clear j n nstates_X nstates_Y options_ode_sim rhs_Y t tmax_epi_rand x toggle_save
    save(['~/Desktop/host_multipath_davies.mat'])
end

function [rhs_X,w_X] = longtime_map(w_X,tags_Ext_NZ,tags_Ext_Y,nX,nY)
    if isempty(nX)
        nX = ones(1,size(tags_Ext_NZ,2));
    end
    if isempty(nY)
        nY = ones(1,size(tags_Ext_Y,2));
    end

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