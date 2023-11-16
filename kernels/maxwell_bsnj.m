function bsnj = maxwell_bsnj(s,n,jj,vts_parallel,Ts_parallel,Ts_perp,wcs,kz,J)
% filename: maxwell_bsnj.m
% Calculate the coefficients of bsnl for the oblique plasma wave model 
% with a bi-Maxweillian distrubution.
% The script created by Bai Wei (email:baiwei12@mail.ustc.edu.cn) 
% was created on Sep 3rd, 2023，
% and revised on Sep 4th, 2023.


% J-pole
[bj,cj] = func_Jpole(J);

bsnj = bj(jj)*cj(jj)*kz*vts_parallel(s)*Ts_perp(s)/Ts_parallel(s) + bj(jj)*n*wcs(s);