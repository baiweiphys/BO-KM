function bz31snj = maxwell_bz31snj(s,n,jj,bsnj,csnj,wps,wcs,lambdas,theta)
% @Description: Calculate the coefficients of bz31snj for the oblique 
% plasma wave model with a bi-Maxwellian distribution.
% @Filename: maxwell_bz31snj.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

params_with_unit;

In = besseli(n,lambdas(s));

coef = wps(s)^2*exp(-1*lambdas(s))/lambdas(s)*In*bsnj(s,n,jj)...
       *(csnj(s,n,jj)/wcs(s)^2 - 2*n/wcs(s) + n^2/csnj(s,n,jj));

bz31snj = -1i*epsilon0*tan(theta)^2*coef;


