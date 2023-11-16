function bx2sn = biKappa_bx2sn(s,n,kappas,wps,lambdas,SnFunc)
% filename: biKappa_bx2sn.m
% Calculate the coefficients of bx2sn for x-component of perpendicular
% plasma waves with a bi-kappa distrubution.
% Generated on October 1st, 2023 by Bai Wei (baiweiphys@gmail.com).
% Modified on Sep 24th, 2023

params_with_unit;

bx2sn = -1*epsilon0*wps(s).^2*n/lambdas(s)^2*SnFunc(s,n);

end
