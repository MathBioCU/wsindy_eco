function [rhs_xy,W,tags_param] = longtime_map(W,lib_state,lib_param,n_state,n_param)
    if isequal(class(lib_state),'double')
        lib_state = library('tags',lib_state);
    end
    if isequal(class(lib_param),'double')
        lib_param = library('tags',lib_param);
    end
    tags_param = cell2mat(lib_param.tags(:));
    tags_state = cell2mat(lib_state.tags(:));
    supp = cellfun(@(W)find(any(W~=0,2)),W(:),'uni',0);
    features_W = cell(length(W),1);
    for i=1:length(W)
        for j=1:length(find(supp{i}))
            f = @(X) X(1)*0;
            for k=1:size(tags_param,1)
                W{i}(supp{i}(j),k) = [W{i}(supp{i}(j),k)/prod(n_state.^lib_state.tags{supp{i}(j)})*n_state(i)]*prod((1./n_param).^tags_param(k,:));
                f = @(X) f(X)+W{i}(supp{i}(j),k)*lib_param.terms{k}.evalterm(X);%prod(X.^tags_NZ(k,:));
            end
            features_W{i}{j} = f;
        end
    end    
    features = cellfun(@(s)lib_state.get_fHandles(s),supp,'uni',0);
    % features = cellfun(@(s)arrayfun(@(i)@(X)prod(X.^tags_state(i,:)),s,'uni',0),supp,'uni',0);
    param_map = @(Y)cellfun(@(f)arrayfun(@(i)f{i}(Y),1:length(f)),features_W,'uni',0);
    % rhs_xy = @(X,Y) rhs_fun_1inp(features,param_map(Y),X);
    rhs_xy = @(y,X) rhs_fun(features,param_map(X),y); % 
    tags_param = cell2mat(lib_param.tags(:));
end

function dx = rhs_fun_1inp(features,params,x)
    dx = zeros(length(params),1);
    for i=1:length(params)
        if ~isempty(features{i})
            dx(i) = dot(cellfun(@(z1) z1(x),features{i}),params{i});
        end
    end
end