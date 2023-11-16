function bx31snl = km_bx31snl(s,n,l,csn,bsl,bsnl,wps,lambdas,wcs,th)
% @Description: Calculate the coefficients of bx31snl for the x-component 
% of oblique plasma waves with a kappa-Maxwellian distribution.
% @Filename: km_bx31snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

params_with_unit;

In = besseli(n,lambdas(s));

bx31snl = -1i*epsilon0*wps(s).^2*tan(th)*exp(-1*lambdas(s))/lambdas(s)...
          *n*In*(bsnl(s,n,l)*(csn(s,n)-n*wcs(s))+bsl(s,l))/wcs(s);
end