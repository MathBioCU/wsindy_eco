function n=peak_diff(X_pred,X_test,tol)
    if size(X_pred,1)<3
        n=0;
    else
        [~,p_pred]=findpeaks(X_pred(:,1));
        [~,p_test]=findpeaks(X_test(:,1));
        L = min(length(p_pred),length(p_test));
        n = find(abs(p_pred(1:L)-p_test(1:L))>tol,1);
        if isempty(n)
            n = L;%size(X_pred,1);
        elseif n==1
            n = L;%size(X_pred,1);
        else
            n = n-1;
        end
    end
end