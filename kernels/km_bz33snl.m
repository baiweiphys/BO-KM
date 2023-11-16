function bz33snl = km_bz33snl(s,n,l,csn,bsl,bsnl,wps,lambdas,wcs,th)
% @Description: Calculate the coefficients of bz33snl for the z-component 
% of oblique plasma waves with a kappa-Maxwellian distribution.
% @Filename: km_bz33snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

params_with_unit;

In = besseli(n,lambdas(s));

bz33snl = -1i*epsilon0*wps(s).^2*tan(th)^2*exp(-1*lambdas(s))/lambdas(s)*In...
          *(2*bsnl(s,n,l)*(csn(s,n)-n*wcs(s)) + bsl(s,l))/wcs(s)^2;
end