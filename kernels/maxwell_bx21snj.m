function bx21snj = maxwell_bx21snj(s,n,jj,bsnj,csnj,wps,lambdas)
% @Description: Calculate the coefficients of bx21snj for the oblique 
% plasma wave model with a bi-Maxwellian distribution.
% @Filename: maxwell_bx21snj.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

params_with_unit;

In = besseli(n,lambdas(s));
dIn = 0.5*besseli(n+1,lambdas(s)) + 0.5*besseli(n-1,lambdas(s));

coef = wps(s)^2*exp(-1*lambdas(s))*n*(In-dIn)*bsnj(s,n,jj)/csnj(s,n,jj);
bx21snj = -1*epsilon0*coef;
