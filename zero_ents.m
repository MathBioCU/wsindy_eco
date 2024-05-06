% b = [1 2 6];
% a = {[1 0 2 3 0 4],[0 5 0 6 0 0 7]};
% anew = zero_ents(a,b);
% 
function a = zero_ents(a,b)
    N  = length(b);
    s = cellfun(@(aa)find(aa~=0),a,'Un',0);
    for i=1:N
        ind = 1;
        L = length(s{ind});
        check = true;
        while check
            if b(i) <= L
                a{ind}(s{ind}(b(i)-L+length(s{ind}))) = 0;
                check = false;
            else
                ind = ind+1;
                L = L + length(s{ind});
            end
        end
    end
end
    