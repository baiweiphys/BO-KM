function bz3 = biKappa_bz3(S,N,kappas,Ts_parallel,Ts_perp,wps,lambdas,SnFunc)
% filename: biKappa_bz3.m
% Calculate the coefficients of bz3 for z-component of perpendicular
% plasma waves with a bi-kappa distrubution.
% Generated on October 1st, 2023 by Bai Wei (baiweiphys@gmail.com).
% Modified on Sep 24th, 2023

params_with_unit;

As = @(s) Ts_perp(s)/Ts_parallel(s) -1.0;
tmp = 0;

Nvector = -N:N;
for s=1:S
    for in=1:length(Nvector)
        n = Nvector(in);
        tmp = tmp + 1i*epsilon0*wps(s).^2/lambdas(s)*As(s)/(As(s)+1.0)*SnFunc(s,n);
    end
end

bz3 = tmp;
              
