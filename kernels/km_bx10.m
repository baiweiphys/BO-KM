function bx10 = km_bx10(wps)
% filename:km_bx10.m
% Calculate the coefficients of bx10 for x-component of oblique 
% plasma wave model with a kappa-Maxweillian distrubution.
% Modified on Sep 24th, 2023

params_with_unit;

bx10 = 1i*epsilon0*sum(wps.^2);
end