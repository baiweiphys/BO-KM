function csnj = maxwell_csnj(s,n,jj,vts_parallel,wcs,us0,kz,J)
% @Description: Calculate the coefficients for the oblique plasma 
% wave model with a bi-Maxwellian distribution.
% @Filename: maxwell_csnj.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

% J-pole
[bj,cj] = func_Jpole(J);

csnj = n*wcs(s) + kz*us0(s) + cj(jj)*kz*vts_parallel(s);
