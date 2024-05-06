%%% evo mod:
b = 0.129; 
lam = 0.21; % 0.21
s = 1.21; 
V = 2.97; 
delta = 1/16; 
m = 27;
a = 0.96; 
w = 0.14; 
phi = 7.4; % 7.4 
gam = 0.3; % 0.3
mu = 0.01; %???
sig = 0; %???
X0s = [0.01 0.001 2];
% x0s = [0.1 0.05 1.57];

%%% simpler mod:
sig = 0;
b = 0; 
s = 0; 
a = 0; 
delta = 1/16; 
mu = 0.01;
m = 0;
lam = 5.5;
V = 0.5;%2.97; *2
phi = 7; 
gam = 0.3;
% X0s = [0.007127252309866   0.018345638426356   1.000000000000000];
X0s = [rand*0.1 rand 1];
% X0s = [0.009903600000741   0.000000050011954 1];
% X0s = [0.007496455808268  0.000014499661094 1]; % 2 before the host peak
% X0s = [0.046992571460427   0.000077174902426 1]; % 2 before pathogen peak

tol_ode = 10^-12;

M = 81;
num_t_epi = 56*2;

fN = @(N,Z,nu,I) lam*exp(sig*randn)*N.*(1-I).*(1-2*a*w*N./(w^2+N^2)).*(1+s*nu.*(1-I).^V);
fZ = @(N,Z,nu,I) phi*N.*I+gam*Z;
fnu = @(N,Z,nu,I) nu.*(1-I).^(b*V).*(1+s*nu*(b*V+1).*(1-I).^(b*V))./(1+s*nu.*(1-I).^(b*V));

N = X0s(1); Z = X0s(2); nu = X0s(3);
Ycell = cell(M-1,1);
t_epi = cell(M-1,1);
yearlength = 56;
sig_tmax = 0;
for j=1:M-1
    tmax_epi_rand = yearlength + 2*(rand-0.5)*sig_tmax; 
    t_epi{j} = linspace(0,tmax_epi_rand,num_t_epi)';
    X0 = [N(end) zeros(1,m) Z(end)];
    [t,y] = get_epi(t_epi{j},X0,nu(end),V,m,delta,mu,tol_ode);
    I = (X0(1)-y(end,1))/X0(1);
    % I = ifrac([x0(1) x0(end)],V,nu(end));
    N(end+1) = fN(X0(1),X0(end),nu(end),I);
    Z(end+1) = fZ(X0(1),X0(end),nu(end),I);
    nu(end+1) = fnu(X0(1),X0(end),nu(end),I);
    % Ycell{j} = y;
    Ycell{j} = y(:,[1 end]);
end
tn = (0:M-1)*yearlength;
figure(1)
semilogy([N;Z]','o-')
legend
X = [N;Z]';

Y = cell2mat(Ycell); 
t = cell2mat(arrayfun(@(i)(i-1)*yearlength+t_epi{i},(1:M-1)','uni',0));
for j=1:2
    subplot(2,1,j)
    plot(tn,X(:,j),'b-o',t,Y(:,j),'r-','linewidth',2)
    legend({'X','Y'})
end

nstates_Y = size(Ycell{1},2);
nstates_X = 2;
E = eye(nstates_Y);
tags_IC_true = [[1 0];[0 1]];
W_IC_true = {[1 zeros(1,nstates_Y-1)],[zeros(1,nstates_Y-1) 1]};
tags_X_true = [[0 0];[-V 0]];
tags_Y_true = [[1+V zeros(1,nstates_Y-2) 1];E(2:end,:)];
if m>0
    W_Y_true = {[[0 -nu(1)];zeros(size(tags_Y_true,1)-1,2)],...
        [[0 nu(1)];[-m*delta 0];zeros(size(tags_Y_true,1)-1,2)]};
    for j=3:nstates_Y-1
        W_Y_true{j} = [zeros(j-2,2);[m*delta 0];[-m*delta 0];zeros(size(tags_Y_true,1)-j,2)];
    end
    W_Y_true{nstates_Y} = [zeros(size(tags_Y_true,1)-2,2);[m*delta 0];[-mu 0]];
else
    W_Y_true = {[[0 -nu(1)];[0 0]],[[0 nu(1)];[-mu 0]]};
end
tags_Ext_X_true = [[0 0];[0 1];[1 0]];
tags_Ext_Y_true = [E(1,:)*0;E(1,:)];
W_X_true = {[[0 lam];[0 0];[0 0]],[[0 -phi];[gam 0];[phi 0]]};

custom_tags_Y = {[1+V zeros(1,nstates_Y-2) 1]}; %<<< decision needs to made
custom_tags_X = {[-V 0]}; %<<< decision needs to made

% linregargs_fun_IC = @(WS){'Aineq',-WS.G{1},'bineq',zeros(size(WS.G{1},1),1)};
% linregargs_fun_Y = @(WS){};
% linregargs_fun_X = @(WS){'Aineq',-WS.G{1},'bineq',zeros(size(WS.G{1},1),1)};

% N = 50; UB = 0.02;
% linregargs_fun_IC = @(WS)enforce_pos_zero(WS);
% linregargs_fun_Y = @(WS)enforce_pos_zero(WS);
% linregargs_fun_X = @(WS)enforce_pos(WS,N,UB);

linregargs_fun_IC = @(WS){};
linregargs_fun_Y = @(WS){};
linregargs_fun_X = @(WS){};

[rhs_IC_true,W_IC,tags_X_IC] = shorttime_map(W_IC_true,library('tags',zeros(1,nstates_Y)),tags_IC_true,ones(1,nstates_Y),ones(1,nstates_X));
rhs_IC_true = @(X) rhs_IC_true(zeros(nstates_Y,1),X(:));
[rhs_Y_true,W_Y,tags_Y_Yeq] = shorttime_map(W_Y_true,tags_Y_true,tags_X_true,ones(1,nstates_Y),ones(1,nstates_X));
[rhs_X_true,W_X] = longtime_map(W_X_true,tags_Ext_X_true,tags_Ext_Y_true,ones(1,nstates_X),ones(1,nstates_Y));

% N=5000;
% 
% subplot(2,1,1)
% plot(tn(1:N)/yearlength,X(1:N,1),'b-o','linewidth',2)
% legend({'N'})
% xlabel('n')
% 
% subplot(2,1,2)
% plot(tn(1:N)/yearlength,X(1:N,2),'b-o','linewidth',2)
% legend({'Z'})
% xlabel('n')
% 
% saveas(gcf,['~/transients_',num2str(N),'.png'])


function out = enforce_pos(WS,N,UB,etc)
   if ~exist('etc','var')
       etc = {};
   end
   N = min(N,size(WS.G{1},1));
   X = lhsdesign(N,WS.nstates)*UB;
   A = [];
   for j=1:WS.numeq
       A = blkdiag(A,WS.lib(j).evalterms(X));
   end
   b = zeros(size(A,1),1);
   out = [{'Aineq',-A,'bineq',b},etc];
end

function out = enforce_pos_zero(WS,etc)
   if ~exist('etc','var')
       etc = {};
   end
   X = zeros(1,WS.nstates);
   A = [];
   for j=1:WS.numeq
       A = blkdiag(A,WS.lib(j).evalterms(X));
   end
   b = zeros(size(A,1),1);
   out = [{'Aineq',-A,'bineq',b},etc];
end

function [t,x] = get_epi(t,x0,nu,V,m,delta,mu,tol_ode)
    options_ode_sim = odeset('RelTol',tol_ode,'AbsTol',tol_ode*ones(1,length(x0)));
    [t,x] = ode45(@(t,x)rhs(x,nu,V,m,delta,mu,x0(1)),t,x0,options_ode_sim);
end

function dx = rhs(X,nu,V,m,delta,mu,S0)
    dx = X*0;
    dx(1) = -nu*X(1).*X(end).*(X(1)/S0).^V;
    if m>0
        dx(2) = nu*X(1).*X(end).*(X(1)/S0).^V - m*delta*X(2);
        for j=3:length(X)-1
            dx(j) = m*delta*X(j-1)-m*delta*X(j);
        end
        dx(end) = m*delta*X(end-1)-mu*X(end);
    else
        dx(2) = nu*X(1).*X(end).*(X(1)/S0).^V - mu*X(end);
    end
end

function I = ifrac(x,C2,nu)
    x1 = x(1);x2 = x(2);
    F = @(I) (1-I).*(1 + (nu*C2*x1)*I + (nu*C2*x2)).^(1/C2) - 1;
    L = (1-(1+nu*C2*x2)/(nu*x1))/(1+C2);
    if L<0
        I0s = [0 1];
    else
        I0s = [L 1];
    end
    I = fzero(F,I0s);
end