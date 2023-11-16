function csn = km_csn(s,n,kappas,vts_parallel,wcs,us0,kz)
% @Description: Calculate the coefficients of csn(s) for the x-component 
% of oblique plasma waves with a kappa-Maxwellian distribution.
% @Filename: km_bx_KM.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15


csn = n*wcs(s) + kz*us0(s) -1i*sqrt(kappas(s))*kz*vts_parallel(s);


