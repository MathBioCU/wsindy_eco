%%% broken script!!! %%%
%%% goal is to estimate T, using the gradient w/r/t T

iter_max = 1000;
alphaT = 50000;

X_Yend = zeros(max(train_inds),nstates_X+nstates_Y);
X_Yend(train_inds,1:nstates_X) = X_train;
subinds = train_inds(diff(train_inds)==1);
polys_X = 0:1;
polys_Y = 0:1;
tags_X_Xeq = get_tags(polys_X,[],nstates_X);
tags_Y_Xeq = get_tags(polys_Y,[],nstates_Y);
lib_X_Xeq = library('tags',tags_X_Xeq);
lib_Y_Xeq = library('tags',tags_Y_Xeq);
lib_Xeq = kron_lib(lib_X_Xeq,lib_Y_Xeq);
E = eye(nstates_X+nstates_Y);

Tends = subinds*0+56;
for iter=1:iter_max
Yend = cell2mat(arrayfun(@(i)interp1(Y_ns{i,2},Y_ns{i,1},Tends(i))./nY,(1:length(subinds))','uni',0));
X_Yend(subinds,nstates_X+1:end) = Yend;
Uobj_X_Yend = wsindy_data(X_Yend,(0:max(train_inds)-1)*yearlength);
tf_X = testfcn(Uobj_X_Yend,'meth','direct','param',1/2,'phifuns','delta','mtmin',0,'subinds',subinds);
WS_Xeq = wsindy_model(Uobj_X_Yend,lib_Xeq,tf_X,'lhsterms', arrayfun(@(i)term('ftag',E(i,:),'linOp',1),(1:nstates_X)','uni',0));
WS_Xeq.cat_Gb('cat','blkdiag');

linregargs_X = {};%linregargs_fun_X(WS_Xeq.G{1},WS_Xeq.b{1}); 
[WS_Xeq,loss_X,lambda_X,w_its,res_X,res_0_X,CovW_X] = WS_opt().MSTLS_WENDy(WS_Xeq,...
    'maxits_wendy',0,'lambdas',0,'ittol',wendy_ittol,'diag_reg',diag_reg,'verbose',wendy_verb,'linregargs',linregargs_X);

L=WS_Xeq.add_L1;
L = L{1}(:,WS_Xeq.dat.dims*nstates_X+1:end);
ff = zeros(WS_Xeq.dat.dims,nstates_Y);
foo = arrayfun(@(i)1./nY'.*rhs_Y(X_Yend(i,nstates_X+1:end).*nY,X_Yend(i,1:nstates_X).*nX),...
    subinds','uni',0);
ff(subinds,:) = cell2mat(foo)';
foo = L*ff(:);
foo = foo.*WS_Xeq.res;
foo = reshape(foo,[],2);
grad_T = -foo*[1;1];
tt=cellfun(@(t)t(end),t_epi);

Tends = Tends-alphaT*grad_T;
disp([Tends tt(subinds)])

end
% G= WS_Xeq.G{1};
% r= WS_Xeq.res;
% grad_W = -G'*r;



