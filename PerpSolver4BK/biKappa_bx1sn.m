function bx1sn = biKappa_bx1sn(s,n,kappas,wps,lambdas,SnFunc)
% filename: biKappa_bx1sn.m
% Calculate the coefficients of bx1sn for x-component of perpendicular
% plasma waves with a bi-kappa distrubution.
% Generated on October 1st, 2023 by Bai Wei (baiweiphys@gmail.com).
% Modified on Sep 24th, 2023

params_with_unit;

bx1sn = 1i*epsilon0*wps(s).^2*n^2/lambdas(s)^2*SnFunc(s,n);

end
