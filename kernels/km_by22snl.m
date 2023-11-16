function by22snl = km_by22snl(s,n,l,bsl,wps,lambdas)
% @Description: Calculate the coefficients of by22snl for the y-component 
% of oblique plasma waves with a kappa-Maxwellian distribution.
% @Filename: coef_by22snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

params_with_unit;

dIn = 0.5*besseli(n+1,lambdas(s)) + 0.5*besseli(n-1,lambdas(s));
In = besseli(n,lambdas(s));

by22snl = -1i*epsilon0*wps(s).^2*exp(-1*lambdas(s))/lambdas(s)*bsl(s,l)...
          *(n^2*In + 2*lambdas(s)^2*(In-dIn));
end