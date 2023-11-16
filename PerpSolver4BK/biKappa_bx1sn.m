function bx1sn = biKappa_bx1sn(s,n,kappas,wps,lambdas,SnFunc)
% @Description: Calculate the coefficients of bx1sn for the x-component of
% perpendicular plasma waves exhibiting a bi-kappa distribution.
% @Filename: biKappa_bx1sn.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2021-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

params_with_unit;

bx1sn = 1i*epsilon0*wps(s).^2*n^2/lambdas(s)^2*SnFunc(s,n);

end
