function bx30 = km_bx30(wps,th)
% @Description: Calculate the coefficients of bx30 for the x-component 
% of oblique plasma waves with a kappa-Maxwellian distribution.
% @Filename: km_bx30.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

params_with_unit;

bx30 = -1i*epsilon0*sum(wps.^2)*tan(th);
end