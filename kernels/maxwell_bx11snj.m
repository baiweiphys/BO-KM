function bx11snj = maxwell_bx11snj(s,n,jj,bsnj,csnj,wps,lambdas)
% @Description: Calculate the coefficients of bx11snj for the oblique 
% plasma wave model with a bi-Maxwellian distribution.
% @Filename: maxwell_bx11snj.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

params_with_unit;

In = besseli(n,lambdas(s));

bx11snj = -1i*epsilon0*wps(s)^2*exp(-1*lambdas(s))/lambdas(s)...
          *n^2*In*bsnj(s,n,jj)/csnj(s,n,jj);


