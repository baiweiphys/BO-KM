function bz31snj = maxwell_bz31snj(s,n,jj,bsnj,csnj,wps,wcs,lambdas,theta)
% filename: maxwell_bz31snj.m
% Calculate the coefficients of bz31snj for the oblique plasma wave model 
% with a bi-Maxweillian distrubution.
% The script created by Bai Wei (email:baiwei12@mail.ustc.edu.cn) 
% was created on Sep 3rd, 2023，
% and revised on Sep 4th, 2023.

params_with_unit;

In = besseli(n,lambdas(s));

coef = wps(s)^2*exp(-1*lambdas(s))/lambdas(s)*In*bsnj(s,n,jj)...
       *(csnj(s,n,jj)/wcs(s)^2 - 2*n/wcs(s) + n^2/csnj(s,n,jj));

bz31snj = -1i*epsilon0*tan(theta)^2*coef;


