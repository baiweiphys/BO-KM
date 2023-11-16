function bx12snl = km_bx12snl(s,n,l,bsl,wps,lambdas)
% @Description: Calculate the coefficients of bx12snl for the x-component 
% of oblique plasma waves with a kappa-Maxwellian distribution.
% @Filename: km_bx12snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

params_with_unit;

In = besseli(n,lambdas(s));

bx12snl = -1i*epsilon0*wps(s).^2*exp(-1*lambdas(s))/lambdas(s)*n^2*In*bsl(s,l);
end