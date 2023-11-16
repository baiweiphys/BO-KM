function bz30 = km_bz30(wps,th)
% @Description: Calculate the coefficients of bz30 for the z-component 
% of oblique plasma waves with a kappa-Maxwellian distribution.
% @Filename: km_bz30.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

params_with_unit;

bz30 = 1i*epsilon0*sum(wps.^2)*tan(th)^2;
end