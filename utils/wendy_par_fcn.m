function [rhs,W,WS,lib_param,W_its,res_WENDy,res_0,CovW] = ...
    wendy_par_fcn(lib_state,lib_param,Uobj,tf,lhs,WENDy_args,linregargs,scale_state,scale_param)

    lib = kron_lib(lib_state,lib_param);
    WS = wsindy_model(Uobj,lib,tf,'lhsterms', lhs);
    WS.cat_Gb('cat','blkdiag');
    WENDy_args = [WENDy_args,{'linregargs',linregargs}];
    [WS,W_its,res_WENDy,res_0,CovW] = WS_opt().wendy(WS,WENDy_args{:});
    W = cellfun(@(w)reshape(w,length(lib_state.terms),[]),WS.reshape_w,'uni',0);
    [rhs,W,~,M] = shorttime_map(W,lib_state,lib_param,scale_state,scale_param);
    CovW = (M.*CovW).*(M');

end