function by31snj = maxwell_by31snj(s,n,jj,bsnj,csnj,wps,wcs,lambdas,theta)
% @Description: Calculate the coefficients of by31snj for the oblique 
% plasma wave model with a bi-Maxwellian distribution.
% @Filename: maxwell_by31snj.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

params_with_unit;

In = besseli(n,lambdas(s));
dIn = 0.5*besseli(n+1,lambdas(s)) + 0.5*besseli(n-1,lambdas(s));

by31snj = epsilon0*tan(theta)*wps(s)^2*exp(-1*lambdas(s))...
          *(In-dIn)*(1-n*wcs(s)/csnj(s,n,jj))*bsnj(s,n,jj)/wcs(s);


