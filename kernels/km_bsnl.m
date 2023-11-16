function bsnl = km_bsnl(s,n,l,kappas,vts_parallel,wcs,kz)
% @Description: Calculate the coefficients of bsnl for the oblique 
% plasma waves with a kappa-Maxwellian distribution.
% @Filename: km_bsnl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

coefbs0 =  @(s,l) 1i/vts_parallel(s)*(kappas(s)-0.5)/sqrt(kappas(s));
coefbs1 =  @(s,l) coefbs0(s,l)*factorial(kappas(s))/factorial(2*kappas(s))...
          *factorial(2*kappas(s)-l-1)/factorial(kappas(s)-l);
coefbs2 = @(s,l) coefbs1(s,l)*(2i*sqrt(kappas(s))*kz*vts_parallel(s)).^l;

bsnl = coefbs2(s,l)*2.0*n*wcs(s)/kz;