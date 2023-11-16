function csnj = maxwell_csnj(s,n,jj,vts_parallel,wcs,us0,kz,J)
% filename: maxwell_csnj.m
% Calculate the coefficients for the oblique plasma wave model 
% with a bi-Maxweillian distrubution.
% The script created by Bai Wei (email:baiwei12@mail.ustc.edu.cn) 
% was created on Sep 3rd, 2023，
% and revised on Oct 16th, 2023.

% J-pole
[bj,cj] = func_Jpole(J);

csnj = n*wcs(s) + kz*us0(s) + cj(jj)*kz*vts_parallel(s);
