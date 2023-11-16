function bz3sn = biKappa_bz3sn(s,n,kappas,Ts_parallel,Ts_perp,wps,lambdas,SnFunc)
% filename: biKappa_bz3sn.m
% Calculate the coefficients of bz3sn for z-component of perpendicular
% plasma waves with a bi-kappa distrubution.
% Generated on October 1st, 2023 by Bai Wei (baiweiphys@gmail.com).
% Modified on Sep 24th, 2023

params_with_unit;

As = Ts_perp/Ts_parallel -1.0;

bz3sn = 1i*epsilon0*wps(s).^2/lambdas(s)/(As(s)+1.0)*SnFunc(s,n);

end