function bx31snj = maxwell_bx31snj(s,n,jj,bsnj,csnj,wps,wcs,lambdas,theta)
% @Description: Calculate the coefficients of bx31snj for the oblique 
% plasma wave model with a bi-Maxwellian distribution.
% @Filename: maxwell_bx31snj.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

params_with_unit;

In = besseli(n,lambdas(s));

coef = wps(s)^2*exp(-1*lambdas(s))/lambdas(s)...
          *n*In*(1-n*wcs(s)/csnj(s,n,jj))*bsnj(s,n,jj)/wcs(s);

bx31snj = -1i*epsilon0*tan(theta)*coef;


