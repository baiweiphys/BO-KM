function bx31snj = maxwell_bx31snj(s,n,jj,bsnj,csnj,wps,wcs,lambdas,theta)
% filename: maxwell_bx31snj.m
% Calculate the coefficients of bx31snj for the oblique plasma wave model 
% with a bi-Maxweillian distrubution.
% The script created by Bai Wei (email:baiwei12@mail.ustc.edu.cn) 
% was created on Sep 3rd, 2023，
% and revised on Sep 4th, 2023.

params_with_unit;

In = besseli(n,lambdas(s));

coef = wps(s)^2*exp(-1*lambdas(s))/lambdas(s)...
          *n*In*(1-n*wcs(s)/csnj(s,n,jj))*bsnj(s,n,jj)/wcs(s);

bx31snj = -1i*epsilon0*tan(theta)*coef;


