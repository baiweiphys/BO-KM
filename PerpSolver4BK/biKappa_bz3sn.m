function bz3sn = biKappa_bz3sn(s,n,kappas,Ts_parallel,Ts_perp,wps,lambdas,SnFunc)
% @Description: Calculate the coefficients of bz3sn for the z-component of
% perpendicular plasma waves exhibiting a bi-kappa distribution.
% @Filename: biKappa_bz3sn.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2021-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

params_with_unit;

As = Ts_perp/Ts_parallel -1.0;

bz3sn = 1i*epsilon0*wps(s).^2/lambdas(s)/(As(s)+1.0)*SnFunc(s,n);

end