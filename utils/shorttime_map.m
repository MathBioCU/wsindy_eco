function [rhs_xy,W,tags_param,M] = shorttime_map(W,lib_state,lib_param,n_state,n_param)
    if isequal(class(lib_state),'double')
        lib_state = library('tags',lib_state);
    end
    if isequal(class(lib_param),'double')
        lib_param = library('tags',lib_param);
    end
    if isempty(n_state)
        n_state = ones(1,lib_state.nstates);
    end
    if isempty(n_param)
        n_param = ones(1,lib_param.nstates);
    end
    
    tags_param = cell2mat(lib_param.tags(:));
    supp = cellfun(@(W)find(any(W~=0,2)),W(:),'uni',0);
    features_W = cell(length(W),1);
    M = [];
    for i=1:length(W)
        for j=1:length(find(supp{i})) % rows / state fcns
            f = @(X) X(1)*0;
            supp_j = find(W{i}(supp{i}(j),:));
            for k=1:length(supp_j) % cols / param fcns
                M = [M; 1/prod(n_state.^lib_state.tags{supp{i}(j)})*n_state(i)*prod((1./n_param).^tags_param(supp_j(k),:))];
                W{i}(supp{i}(j),supp_j(k)) = W{i}(supp{i}(j),supp_j(k))*M(end);
                f = @(X) f(X)+W{i}(supp{i}(j),supp_j(k))*lib_param.terms{supp_j(k)}.evalterm(X);
            end
            features_W{i}{j} = f;
        end
    end
    
    features = cellfun(@(s)lib_state.get_fHandles(s),supp,'uni',0);
    param_map = @(X)cellfun(@(f)arrayfun(@(i)f{i}(X),1:length(f)),features_W,'uni',0);
    rhs_xy = @(y,X) rhs_fun(features,param_map(X),y); % 
    tags_param = cell2mat(lib_param.tags(:));
end
