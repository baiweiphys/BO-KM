function bz34snl = km_bz34snl(s,n,l,bsnl,wps,lambdas,wcs,th)
% filename: km_bz34snl.m
% Calculate the coefficients of bz34snl for z-component of oblique 
% plasma wave model with a kappa-Maxweillian distrubution.
% Modified on Sep 24th, 2023

params_with_unit;

In = besseli(n,lambdas(s));

bz34snl = -1i*epsilon0*wps(s).^2*tan(th)^2*exp(-1*lambdas(s))/lambdas(s)*In...
          *bsnl(s,n,l)/wcs(s)^2;
end