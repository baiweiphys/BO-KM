function bz30 = km_bz30(wps,th)
% filename: km_bz30.m
% Calculate the coefficients of bz30 for z-component of oblique 
% plasma wave model with a kappa-Maxweillian distrubution.
% Modified on Sep 24th, 2023

params_with_unit;

bz30 = 1i*epsilon0*sum(wps.^2)*tan(th)^2;
end