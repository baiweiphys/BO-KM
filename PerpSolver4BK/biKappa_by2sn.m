function by2sn = biKappa_by2sn(s,n,kappas,wps,lambdas,SnFunc)
% filename: biKappa_by2sn.m
% Calculate the coefficients of by2sn for y-component of perpendicular
% plasma waves with a bi-kappa distrubution.
% Generated on October 1st, 2023 by Bai Wei (baiweiphys@gmail.com).
% Modified on Sep 24th, 2023

params_with_unit;

by2sn = 1i*epsilon0*wps(s).^2/lambdas(s)^2*SnFunc(s,n);

end