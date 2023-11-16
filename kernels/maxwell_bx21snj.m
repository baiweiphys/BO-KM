function bx21snj = maxwell_bx21snj(s,n,jj,bsnj,csnj,wps,lambdas)
% filename: maxwell_bx21snj.m
% Calculate the coefficients of bx21snj for the oblique plasma wave model 
% with a bi-Maxweillian distrubution.
% The script created by Bai Wei (email:baiwei12@mail.ustc.edu.cn) 
% was created on Sep 3rd, 2023，
% and revised on Sep 4th, 2023.

params_with_unit;

In = besseli(n,lambdas(s));
dIn = 0.5*besseli(n+1,lambdas(s)) + 0.5*besseli(n-1,lambdas(s));

coef = wps(s)^2*exp(-1*lambdas(s))*n*(In-dIn)*bsnj(s,n,jj)/csnj(s,n,jj);
bx21snj = -1*epsilon0*coef;
