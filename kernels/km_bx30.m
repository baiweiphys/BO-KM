function bx30 = km_bx30(wps,th)
% filename: km_bx30.m
% Calculate the coefficients of bx30 for x-component of oblique 
% plasma wave model with a kappa-Maxweillian distrubution.
% Modified on Sep 24th, 2023

params_with_unit;

bx30 = -1i*epsilon0*sum(wps.^2)*tan(th);
end