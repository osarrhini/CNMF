function [w,h] = bodvarsson_wh_scaling( w, h )
%bodvarsson_wh_scaling Using Bodvarsson scaling of W for PET data
%   w and h are the current Non negative matrix
% **This code comes with no guarantee or warranty of any kind.**
% ------------------------------------------------------------------------
    [n,r] = size(w);
    m = size(h,2);
    mm = (w'*w);
    inv_mm = mm\eye(r);
    u=inv_mm * (sum(w,1))';
    w = w.*(ones(n,1)*u');
    inv_u = 1./u;
    h = h.*(inv_u * ones(1,m));
end