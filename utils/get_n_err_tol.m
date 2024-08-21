function n_err_tol = get_n_err_tol(X_pred,X_test,err_tol)
    n = min(size(X_pred,1),size(X_test,1));
    if n>0
        cumerr = arrayfun(@(i)norm(vecnorm(X_pred(1:i,:)-X_test(1:i,:),2,2))/norm(vecnorm(X_test(1:i,:),2,2)),(1:n)');
        n_err_tol = find(cumerr>err_tol,1);
        if isempty(n_err_tol)
            n_err_tol = n-1;
        else
            n_err_tol = n_err_tol-1;
        end
    else
        n_err_tol = 0;
    end
end