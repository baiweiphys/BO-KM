function bx32snl = km_bx32snl(s,n,l,csn,bsl,wps,lambdas,wcs,th)
% @Description: Calculate the coefficients of bx32snl for the x-component 
% of oblique plasma waves with a kappa-Maxwellian distribution.
% @Filename: km_bx32snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

params_with_unit;

In = besseli(n,lambdas(s));

bx32snl = -1i*epsilon0*wps(s).^2*tan(th)*exp(-1*lambdas(s))/lambdas(s)...
          *n*In*bsl(s,l)*(csn(s,n)-n*wcs(s))/wcs(s);
end