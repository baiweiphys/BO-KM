function bsnj = maxwell_bsnj(s,n,jj,vts_parallel,Ts_parallel,Ts_perp,wcs,kz,J)
% @Description: Calculate the coefficients of bsnl for the oblique 
% plasma wave model with a bi-Maxwellian distribution.
% @Filename: maxwell_bsnj.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15


% J-pole
[bj,cj] = func_Jpole(J);

bsnj = bj(jj)*cj(jj)*kz*vts_parallel(s)*Ts_perp(s)/Ts_parallel(s) + bj(jj)*n*wcs(s);