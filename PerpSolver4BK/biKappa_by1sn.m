function by1sn = biKappa_by1sn(s,n,kappas,wps,lambdas,SnFunc)
% filename: biKappa_by1sn.m
% Calculate the coefficients of by1sn for y-component of perpendicular
% plasma waves with a bi-kappa distrubution.
% Generated on October 1st, 2023 by Bai Wei (baiweiphys@gmail.com).
% Modified on Sep 24th, 2023

params_with_unit;

by1sn = epsilon0*wps(s).^2*n/lambdas(s)^2*SnFunc(s,n);

end