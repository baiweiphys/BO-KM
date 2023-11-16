function bz21snj = maxwell_bz21snj(s,n,jj,bsnj,csnj,wps,wcs,lambdas,theta)
% filename: maxwell_bz21snj.m
% Calculate the coefficients of bz21snj for the oblique plasma wave model 
% with a bi-Maxweillian distrubution.
% The script created by Bai Wei (email:baiwei12@mail.ustc.edu.cn) 
% was created on Sep 3rd, 2023，
% and revised on Sep 4th, 2023.

params_with_unit;

In = besseli(n,lambdas(s));
dIn = 0.5*besseli(n+1,lambdas(s)) + 0.5*besseli(n-1,lambdas(s));

coef = wps(s)^2*exp(-1*lambdas(s))*(In-dIn)...
       *(1-n*wcs(s)/csnj(s,n,jj))*bsnj(s,n,jj)/wcs(s);

bz21snj = -1*epsilon0*tan(theta)*coef;


