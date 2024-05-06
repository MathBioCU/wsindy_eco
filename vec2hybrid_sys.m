function [W_IC,W_Y,W_X,rhs_IC,rhs_Y,rhs_X] = vec2hybrid_sys(w_vec,ss,...
    tags_IC_X,tags_Y_Yeq,tags_X_Yeq,tags_X_Xeq,tags_Y_Xeq,...
    nstates_X,nstates_Y)
    W_IC_z = arrayfun(@(i)zeros(1,size(tags_IC_X,1)),1:nstates_Y,'Un',0);
    W_Y_z = arrayfun(@(i)zeros(size(tags_Y_Yeq,1),size(tags_X_Yeq,1)),1:nstates_Y,'Un',0);
    W_X_z = arrayfun(@(i)zeros(size(tags_X_Xeq,1),size(tags_Y_Xeq,1)),1:nstates_X,'Un',0);
    W_new = [W_IC_z,W_Y_z,W_X_z];
    
    ind = 0;
    for i=1:length(W_new)
        W_new{i}(ss{i}) = w_vec(ind+(1:length(ss{i})));
        ind = ind+length(ss{i});
    end
    W_IC = W_new(1:2);
    W_Y = W_new(3:4);
    W_X = W_new(5:6);
    
    rhs_IC = shorttime_map(W_IC,zeros(1,nstates_Y),tags_IC_X,ones(1,nstates_Y),ones(1,nstates_X));
    rhs_IC = @(X) rhs_IC(zeros(nstates_Y,1),X(:));
    rhs_Y= shorttime_map(W_Y,tags_Y_Yeq,tags_X_Yeq,ones(1,nstates_Y),ones(1,nstates_X));
    rhs_X = longtime_map(W_X,tags_X_Xeq,tags_Y_Xeq,ones(1,nstates_X),ones(1,nstates_Y));
end