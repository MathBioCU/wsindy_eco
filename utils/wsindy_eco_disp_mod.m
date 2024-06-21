%%% W_X
s_X = get_str(W_X,lib_Y_Xeq,lib_X_Xeq);
s_Y = get_str(W_Y,lib_X_Yeq,lib_Y_Yeq);
s_IC = get_str(W_IC,lib_X_IC,lib_Y_IC);

function s = get_str(W_X,lib_Y_Xeq,lib_X_Xeq)
    s = cell(length(W_X),1);
    for i=1:length(W_X)
        s{i}='';
        for j=1:size(W_X{i},1)
            for k=1:size(W_X{i},2)
                if W_X{i}(j,k)~=0
                    s{i} = [s{i},['(',num2str(W_X{i}(j,k)),')'],lib_Y_Xeq.terms{k}.get_str,lib_X_Xeq.terms{j}.get_str,'+'];
                end
            end
        end
        s{i} = s{i}(1:end-1);
    end
end