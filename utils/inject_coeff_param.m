function wnew = inject_coeff_param(w,tags_old,tags_p_old,tags_new,tags_p_new)
    wnew = cell(size(w));
    for i=1:length(w)
        wnew{i} = zeros(size(tags_new,1),size(tags_p_new,1));
        for j=1:size(tags_old,1)
            for k=1:size(tags_p_old,1)
                try
                    wnew{i}(ismember(tags_new,tags_old(j,:),'rows'),ismember(tags_p_new,tags_p_old(k,:),'rows')) = w{i}(j,k);
                catch
                    disp('term not included in library')
                    continue
                end
            end
        end
    end
end