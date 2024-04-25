%%% solve the adjoint equation in weak form:
%%% cost fcn g = int_0^T S(y(w)) dt, with S(y) = 0.5||y-y_dat||^2
%%% adjoint eq is Psi' = -dF/dy^T*Psi + dS/dy
%%% this equation can be solved in weak form: <-Phi'+Phi*dF/dy,Psi> = <Phi,dS/dy>

%%% get RHS term <Phi,dS/dy>
Uobj_Ynew = cellfun(@(x,s)wsindy_data(x,train_time),Y_new);
Uobj_tot_new = arrayfun(@(i)...
    wsindy_data([[Uobj_Ynew(i).Uobs{:}] repmat(X_train(X_in(i),:),Uobj_Ynew(i).dims,1)],train_time),(1:length(Uobj_Ynew))');
xRHS = WS_Yeq.get_adjoint_rhs(Uobj_tot_new);

%%% get adjoint operator s.t. L*Psi = <Phi,dS/dy>
L = WS_Yeq.get_adjoint_L;

%%% regularize L
gamma = 10^-2;
% regmat = speye(size(L,2));
regmat = spdiags(ones(size(L ,2),1)*[-1 1],0:1,size(L,2)-1,size(L,2));
Lreg = [L;regmat*gamma*norm(L)/norm(full(regmat))];
xRHSreg = [xRHS;zeros(size(regmat,1),1)];

%%% remove end pts (Psi(T) = 0)
Lsub = Lreg*kron(speye(nstates_Y*WS_Yeq.ntraj),[speye(length(train_time)-1);zeros(1,length(train_time)-1)]);

%%% get Psi
Psisub = Lsub\xRHSreg;
Psisub = reshape(Psisub,length(train_time)-1,[]);
Psi = [];
for j=1:nstates_Y*WS_Yeq.ntraj
    Psi = [Psi;Psisub(:,j);0];
end


ind = 5;
Psidat = Psi(nstates_Y*length(train_time)*(ind-1)+1:nstates_Y*length(train_time)*ind);
Psidat = reshape(Psidat,[],nstates_Y);
plot(Psidat,'o-')

S_Y = WS_Yeq.get_supp;
w_Y_sparse = cellfun(@(w,s)w(s),WS_Yeq.reshape_w,S_Y,'uni',0);
partialwf = arrayfun(@(i)arrayfun(@(j)...
    reshape(get_delYdot(Y_new{i}(j,:),X_train(X_in(i),:),WS_Yeq,0,S_Y,w_Y_sparse),[],nstates_Y)',(1:length(train_time))','uni',0)...
    ,(1:length(Y_train))','uni',0);
partialwf = cellfun(@(pw)cell2mat(arrayfun(@(i)cell2mat(cellfun(@(p) p(i,:),pw,'uni',0)),(1:nstates_Y)','uni',0)),partialwf,'uni',0);
partialwf = cell2mat(partialwf);
delG = -Psi'*partialwf*mean(diff(train_time))

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