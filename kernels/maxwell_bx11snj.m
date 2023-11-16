function bx11snj = maxwell_bx11snj(s,n,jj,bsnj,csnj,wps,lambdas)
% filename: maxwell_bx11snj.m
% Calculate the coefficients of bx11snj for the oblique plasma wave model 
% with a bi-Maxweillian distrubution.
% The script created by Bai Wei (email:baiwei12@mail.ustc.edu.cn) 
% was created on Sep 3rd, 2023，
% and revised on Sep 4th, 2023.

params_with_unit;

In = besseli(n,lambdas(s));

bx11snj = -1i*epsilon0*wps(s)^2*exp(-1*lambdas(s))/lambdas(s)...
          *n^2*In*bsnj(s,n,jj)/csnj(s,n,jj);


