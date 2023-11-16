function bx33snl = km_bx33snl(s,n,l,bsnl,wps,lambdas,wcs,th)
% filename: km_bx33snl.m
% Calculate the coefficients of bx33snl for x-component of oblique 
% plasma wave model with a kappa-Maxweillian distrubution.
% Modified on Sep 24th, 2023

params_with_unit;

In = besseli(n,lambdas(s));

bx33snl = -1i*epsilon0*wps(s).^2*tan(th)*exp(-1*lambdas(s))/lambdas(s)...
          *n*In*bsnl(s,n,l)/wcs(s);
end