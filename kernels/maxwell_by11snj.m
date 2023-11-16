function by11snj = maxwell_by11snj(s,n,j,bsnj,csnj,wps,lambdas)
% @Description: Calculate the coefficients of by11snj for the oblique 
% plasma wave model with a bi-Maxwellian distribution.
% @Filename: maxwell_by11snj.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

params_with_unit;

In = besseli(n,lambdas(s));
dIn = 0.5*besseli(n+1,lambdas(s)) + 0.5*besseli(n-1,lambdas(s));

by11snj = epsilon0*wps(s)^2*exp(-1*lambdas(s))*n*(In-dIn)*bsnj(s,n,j)/csnj(s,n,j);
