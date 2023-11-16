function bz34snl = km_bz34snl(s,n,l,bsnl,wps,lambdas,wcs,th)
% @Description: Calculate the coefficients of bz34snl for the z-component 
% of oblique plasma waves with a kappa-Maxwellian distribution.
% @Filename: km_bz34snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

params_with_unit;

In = besseli(n,lambdas(s));

bz34snl = -1i*epsilon0*wps(s).^2*tan(th)^2*exp(-1*lambdas(s))/lambdas(s)*In...
          *bsnl(s,n,l)/wcs(s)^2;
end