function by21snl = km_by21snl(s,n,l,bsnl,wps,lambdas)
% @Description: Calculate the coefficients of by21snl for the y-component 
% of oblique plasma waves with a kappa-Maxwellian distribution.
% @Filename: km_by12snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

params_with_unit;

dIn = 0.5*besseli(n+1,lambdas(s)) + 0.5*besseli(n-1,lambdas(s));
In = besseli(n,lambdas(s));

by21snl = -1i*epsilon0*wps(s).^2*exp(-1*lambdas(s))/lambdas(s)*bsnl(s,n,l)...
          *(n^2*In + 2*lambdas(s)^2*(In-dIn));
end