addpath(genpath('../utils'))
addpath(genpath('../wsindy_obj_base'));

%% feedback FitzHugh-Nagumo

%discrete : N,Z
%continuos: S,P

nstates_Y = 2;
nstates_X = 2;

x0 = [0.007127252309866   0.018345638426356];

%%% evo mod:
w = 0.14; 
sig = 0;
b = 0; 
s = 0; 
a = 0; 
delta = 1/16; 
mu = 0.01;
m = 0;
lam = 5.5;
V = 0.5;
phi = 7; 
gam = 0.3;
nu = 1;

yearlength = 56;
sig_tmax = 0;

tags_IC_true = [[1 0];[0 1]];
W_IC_true = {[1 0],[0 1]};

tags_X_true = [0 0;-V 0];
tags_Y_true = [1+V 1;0 1];
W_Y_true = {[[0 -nu(1)];[0 0]],[[0 nu(1)];[-mu 0]]};

tags_Ext_X_true = [[0 0];[0 1];[1 0]];
tags_Ext_Y_true = [0 0;1 0];
W_X_true = {[[0 lam];[0 0];[0 0]],[[0 -phi];[gam 0];[phi 0]]};

custom_tags_Y = {[1+V zeros(1,nstates_Y-2) 1]}; %<<< decision needs to made
custom_tags_X = {[-V 0]}; %<<< decision needs to made
linregargs_fun_IC = @(WS){};
linregargs_fun_Y = @(WS){};
linregargs_fun_X = @(WS){};

%% sim

toggle_save = 0;
num_gen = 81;
num_t_epi = yearlength*2;
tol_ode = 10^-12;

tic;
[rhs_IC_true,~,~] = shorttime_map(W_IC_true,zeros(1,nstates_Y),tags_IC_true,[],[]);
rhs_IC_true = @(X) rhs_IC_true(zeros(nstates_Y,1),X(:));
[rhs_Y_true,~,~] = shorttime_map(W_Y_true,tags_Y_true,tags_X_true,[],[]);
[rhs_X_true,~,~] = shorttime_map(W_X_true,tags_Ext_X_true,tags_Ext_Y_true,[],[]);

Ycell = cell(num_gen-1,1);
t_epi = cell(num_gen-1,1);
X = zeros(num_gen,nstates_X);
X(1,:) = x0;
options_ode_sim = odeset('RelTol',tol_ode,'AbsTol',tol_ode*ones(1,nstates_Y));
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
toc

if toggle_save==1
    save(['~/Desktop/forced_FHN.mat'])
end

toc

%%% view
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