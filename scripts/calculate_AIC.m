%% Calculate AIC
dr = '../data/';
load([dr,'UQ_plots_corrected_red_model.mat'],'WS_IC','WS_Xeq','WS_Yeq')

r = 0;num_param=0;diag_reg = 10^-1;
[G_IC,b_IC]=WS_IC.apply_cov(WS_IC.G{1},WS_IC.b{1},diag_reg);
C = WS_IC.cov;
r = r + norm(G_IC*WS_IC.weights-b_IC)^2; %+ log(det((1-diag_reg)*C + diag_reg*speye(size(C,1))));
num_param = num_param + length(find(WS_IC.weights));

[G_Y,b_Y]=WS_Yeq.apply_cov(WS_Yeq.G{1},WS_Yeq.b{1},diag_reg);
C = WS_Yeq.cov;
r = r + norm(G_Y*WS_Yeq.weights-b_Y)^2; %+ log(det((1-diag_reg)*C + diag_reg*speye(size(C,1))));
num_param = num_param + length(find(WS_Yeq.weights));

[G_X,b_X]=WS_Xeq.apply_cov(WS_Xeq.G{1},WS_Xeq.b{1},diag_reg);
C = WS_Xeq.cov;
r = r + norm(G_X*WS_Xeq.weights-b_X)^2; %+ log(det((1-diag_reg)*C + diag_reg*speye(size(C,1))));
num_param = num_param + length(find(WS_Xeq.weights));

AIC = 2*num_param+r;
disp([r 2*num_param AIC])
