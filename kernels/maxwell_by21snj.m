function by21snj = maxwell_by21snj(s,n,jj,bsnj,csnj,wps,lambdas)
% filename: maxwell_by21snj.m
% Calculate the coefficients of by21snj for the oblique plasma wave model 
% with a bi-Maxweillian distrubution.
% The script created by Bai Wei (email:baiwei12@mail.ustc.edu.cn) 
% was created on Sep 3rd, 2023，
% and revised on Sep 4th, 2023.

params_with_unit;

In = besseli(n,lambdas(s));
dIn = 0.5*besseli(n+1,lambdas(s)) + 0.5*besseli(n-1,lambdas(s));

by21snj = -1i*epsilon0*wps(s)^2*exp(-1*lambdas(s))/lambdas(s)...
          *(n^2*In + 2*lambdas(s)^2*(In-dIn))*bsnj(s,n,jj)/csnj(s,n,jj);

