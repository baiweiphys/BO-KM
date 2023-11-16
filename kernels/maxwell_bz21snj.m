function bz21snj = maxwell_bz21snj(s,n,jj,bsnj,csnj,wps,wcs,lambdas,theta)
% @Description: Calculate the coefficients of bz21snj for the oblique 
% plasma wave model with a bi-Maxwellian distribution.
% @Filename: maxwell_bz21snj.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

params_with_unit;

In = besseli(n,lambdas(s));
dIn = 0.5*besseli(n+1,lambdas(s)) + 0.5*besseli(n-1,lambdas(s));

coef = wps(s)^2*exp(-1*lambdas(s))*(In-dIn)...
       *(1-n*wcs(s)/csnj(s,n,jj))*bsnj(s,n,jj)/wcs(s);

bz21snj = -1*epsilon0*tan(theta)*coef;


