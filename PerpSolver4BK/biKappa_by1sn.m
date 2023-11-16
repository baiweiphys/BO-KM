function by1sn = biKappa_by1sn(s,n,kappas,wps,lambdas,SnFunc)
% @Description: Calculate the coefficients of by1sn for the y-component of
% perpendicular plasma waves exhibiting a bi-kappa distribution.
% @Filename: biKappa_by1sn.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2021-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

params_with_unit;

by1sn = epsilon0*wps(s).^2*n/lambdas(s)^2*SnFunc(s,n);

end