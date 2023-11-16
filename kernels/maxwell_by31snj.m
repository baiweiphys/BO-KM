function by31snj = maxwell_by31snj(s,n,jj,bsnj,csnj,wps,wcs,lambdas,theta)
% filename: maxwell_by31snj.m
% Calculate the coefficients of by31snj for the oblique plasma wave model 
% with a bi-Maxweillian distrubution.
% The script created by Bai Wei (email:baiwei12@mail.ustc.edu.cn) 
% was created on Sep 3rd, 2023，
% and revised on Sep 4th, 2023.

params_with_unit;

In = besseli(n,lambdas(s));
dIn = 0.5*besseli(n+1,lambdas(s)) + 0.5*besseli(n-1,lambdas(s));

by31snj = epsilon0*tan(theta)*wps(s)^2*exp(-1*lambdas(s))...
          *(In-dIn)*(1-n*wcs(s)/csnj(s,n,jj))*bsnj(s,n,jj)/wcs(s);


