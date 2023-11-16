function by31snl =km_by31snl(s,n,l,csn,bsl,bsnl,wps,lambdas,wcs,th)
% @Description: Calculate the coefficients of by31snl for the y-component 
% of oblique plasma waves with a kappa-Maxwellian distribution.
% @Filename: km_by31snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

params_with_unit;

dIn = 0.5*besseli(n+1,lambdas(s)) + 0.5*besseli(n-1,lambdas(s));
In = besseli(n,lambdas(s));

by31snl = epsilon0*wps(s).^2*tan(th)*exp(-1*lambdas(s))*(In-dIn)...
          *(bsnl(s,n,l)*(csn(s,n)-n*wcs(s)) + bsl(s,l))/wcs(s);
end