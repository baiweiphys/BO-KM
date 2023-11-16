function by32snl = km_by32snl(s,n,l,csn,bsl,wps,lambdas,wcs,th)
% @Description: Calculate the coefficients of by32snl for the y-component 
% of oblique plasma waves with a kappa-Maxwellian distribution.
% @Filename: km_by32snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

params_with_unit;

dIn = 0.5*besseli(n+1,lambdas(s)) + 0.5*besseli(n-1,lambdas(s));
In = besseli(n,lambdas(s));

by32snl = epsilon0*wps(s).^2*tan(th)*exp(-1*lambdas(s))*(In-dIn)...
          *bsl(s,l)*(csn(s,n)-n*wcs(s))/wcs(s);
end