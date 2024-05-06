
function [value, isterminal, direction] = myEvent(T, Y, thresh, toggle_zero_crossing)
    if toggle_zero_crossing
        value      = or(norm(Y) >= thresh, any(Y==0));
    else
        value      = norm(Y) >= thresh;
    end
    isterminal = 1;   % Stop the integration
    direction  = 0;
end