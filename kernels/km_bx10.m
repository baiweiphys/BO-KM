function bx10 = km_bx10(wps)
% @Description: Calculate the coefficients of bx10 for the x-component 
% of oblique plasma waves with a kappa-Maxwellian distribution.
% @Filename: km_bx10.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

params_with_unit;

bx10 = 1i*epsilon0*sum(wps.^2);
end