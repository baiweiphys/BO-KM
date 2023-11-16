function bz3 = biKappa_bz3(S,N,kappas,Ts_parallel,Ts_perp,wps,lambdas,SnFunc)
% @Description: Calculate the coefficients of bz3 for the z-component of
% perpendicular plasma waves exhibiting a bi-kappa distribution.
% @Filename: biKappa_bz3.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2021-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

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
              
