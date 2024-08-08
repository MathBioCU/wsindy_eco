addpath(genpath('../utils'))
addpath(genpath('../wsindy_obj_base'));

%% two-pathogen system parameters
%discrete : N,Z1,Z2
%continuos: S,P1,P2,nu1,nu2,I1

nstates_Y = 6;
nstates_X = 3;

yearlength = 8;
rho = -0.5;
c1 = 2.06;
c2 = 0.68;
mu1 = 0.99;
mu2 = 0.32;
lam = 5;
phi1 = 2;
phi2 = 2;
nu1 = 5.0;
nu2 = 0.5;

S0 = mu1/nu1;
P0 = log(lam)/nu1/yearlength;
x0 = [S0 P0 10.^-4];

tags_IC_true = [eye(3);zeros(1,3)];
W_IC_true = {[1 0 0 0],[0 1 0 0],[0 0 1 0],[0 0 0 nu1],[0 0 0 nu2],[0 0 0 0]};

tags_X_true = [zeros(1,3)]; 
tags_Y_true = [[1 1 0 1 0 0];[1 0 1 0 1 0];... 
                [0 1 0 2 0 0];[0 0 1 1 1 0];...
                [0 0 1 0 2 0];[0 1 0 1 1 0];...
                [0 1 0 0 0 0];[0 0 1 0 0 0]];
W_Y_true = {[[-1 -1 zeros(1,6)]'],...
            [[1 zeros(1,5) -mu1 0]'],...
            [[0 1 zeros(1,5) -mu2]'],...
            [[0 0 -c1^2 -rho*c1*c2 0 0 0 0]'],...
            [[0 0 0 0 -c2^2 -rho*c1*c2 0 0]'],...
            [1;zeros(7,1)],...
            };

tags_Ext_X_true = [0 0 0;...
                    1 0 0];
tags_Ext_Y_true = [zeros(1,6);...
                   1 zeros(1,5);...
                   zeros(1,5) 1;...
                   ];
W_X_true = {[0 lam 0;0 0 0],...
            [0 0 phi1;zeros(1,3)],...
            [0 -phi2 -phi2;phi2 0 0],...
            };

custom_tags_Y = [];
custom_tags_X = [];
linregargs_fun_IC = @(G,b){};%{'Aineq',-G,'bineq',zeros(size(G,1),1)};
linregargs_fun_Y = @(G,b){};
linregargs_fun_X = @(G,b){};

%% sim 
num_gen = 10;
num_t_epi = 224;
tol_ode = 10^-12;
sig_tmax = 0;
toggle_save = 0;

tic
[rhs_IC_true,~,~] = shorttime_map(W_IC_true,zeros(1,nstates_Y),tags_IC_true,[],[]);
rhs_IC_true = @(X) rhs_IC_true(zeros(nstates_Y,1),X(:));
[rhs_Y_true,~,~] = shorttime_map(W_Y_true,tags_Y_true,tags_X_true,[],[]);
[rhs_X_true,~,~] = shorttime_map(W_X_true,tags_Ext_X_true,tags_Ext_Y_true,[],[]);

Ycell = cell(num_gen-1,1);
t_epi = cell(num_gen-1,1);
X = zeros(num_gen,nstates_X);
X(1,:) = x0;
options_ode_sim = odeset('RelTol',tol_ode,'AbsTol',tol_ode*ones(1,nstates_Y));
n=1;
while n<num_gen-1
    disp(n)
    if or(and(X(n,2)>=0,X(n,3)>0),and(X(n,3)>=0,X(n,2)>0))
        tmax_epi_rand = yearlength + 2*(rand-0.5)*sig_tmax;
        t_epi{n} = linspace(0,tmax_epi_rand,num_t_epi)';
        rhs_Y = @(y)rhs_Y_true(y,X(n,:));
        [t_epi{n},x] = ode15s(@(t,x)rhs_Y(x),t_epi{n},rhs_IC_true(X(n,:)),options_ode_sim);
        Ycell{n} = x;
        new_X = rhs_X_true(X(n,:),x(end,:));
        new_X(new_X<tol_ode) = 0;
        X(n+1,:) = new_X;
        n=n+1;
    else
        break
    end
end
if n~=num_gen
    num_gen=n;
    Ycell = Ycell(1:num_gen-1);
    t_epi = t_epi(1:num_gen-1);
    X = X(1:num_gen,:);
end
toc

%%% view 
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

%%% save
if toggle_save==1
    clear j n nstates_X nstates_Y options_ode_sim rhs_Y t tmax_epi_rand x toggle_save
    save(['~/Desktop/host_multipath_6-3_d_steady_state_1path.mat'])
end
