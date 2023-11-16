function by2sn = biKappa_by2sn(s,n,kappas,wps,lambdas,SnFunc)
% @Description: Calculate the coefficients of by2sn for the y-component of
% perpendicular plasma waves exhibiting a bi-kappa distribution.
% @Filename: biKappa_by2sn.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2021-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

params_with_unit;

by2sn = 1i*epsilon0*wps(s).^2/lambdas(s)^2*SnFunc(s,n);

end