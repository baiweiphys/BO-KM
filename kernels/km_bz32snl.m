function bz32snl = km_bz32snl(s,n,l,csn,bsl,wps,lambdas,wcs,th)
% @Description: Calculate the coefficients of bz32snl for the z-component 
% of oblique plasma waves with a kappa-Maxwellian distribution.
% @Filename: km_bz32snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

params_with_unit;

In = besseli(n,lambdas(s));

bz32snl = -1i*epsilon0*wps(s).^2*tan(th)^2*exp(-1*lambdas(s))/lambdas(s)*In...
          *bsl(s,l)*(csn(s,n)-n*wcs(s))^2/wcs(s)^2;
end