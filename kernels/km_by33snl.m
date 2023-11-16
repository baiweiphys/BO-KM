function by33snl = km_by33snl(s,n,l,bsnl,wps,lambdas,wcs,th)
% @Description: Calculate the coefficients of by33snl for the y-component 
% of oblique plasma waves with a kappa-Maxwellian distribution.
% @Filename: km_by32snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

params_with_unit;

% dIn = besseli(n+1,lambdas(s)) + n*besseli(n,lambdas(s))/lambdas(s);
dIn = 0.5*besseli(n+1,lambdas(s)) + 0.5*besseli(n-1,lambdas(s));
In = besseli(n,lambdas(s));

by33snl = epsilon0*wps(s).^2*tan(th)*exp(-1*lambdas(s))*(In-dIn)...
          *bsnl(s,n,l)/wcs(s);
end