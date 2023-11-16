function bz33snl = km_bz33snl(s,n,l,csn,bsl,bsnl,wps,lambdas,wcs,th)
% filename: km_bz33snl.m
% Calculate the coefficients of bz33snl for z-component of oblique 
% plasma wave model with a kappa-Maxweillian distrubution.
% Modified on Sep 24th, 2023

params_with_unit;

In = besseli(n,lambdas(s));

bz33snl = -1i*epsilon0*wps(s).^2*tan(th)^2*exp(-1*lambdas(s))/lambdas(s)*In...
          *(2*bsnl(s,n,l)*(csn(s,n)-n*wcs(s)) + bsl(s,l))/wcs(s)^2;
end