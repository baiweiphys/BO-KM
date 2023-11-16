function by11snj = maxwell_by11snj(s,n,j,bsnj,csnj,wps,lambdas)
% filename: maxwell_by11snj.m
% Calculate the coefficients of by11snj for the oblique plasma wave model 
% with a bi-Maxweillian distrubution.
% The script created by Bai Wei (email:baiwei12@mail.ustc.edu.cn) 
% was created on Sep 3rd, 2023，
% and revised on Sep 4th, 2023.

params_with_unit;

In = besseli(n,lambdas(s));
dIn = 0.5*besseli(n+1,lambdas(s)) + 0.5*besseli(n-1,lambdas(s));

by11snj = epsilon0*wps(s)^2*exp(-1*lambdas(s))*n*(In-dIn)*bsnj(s,n,j)/csnj(s,n,j);
