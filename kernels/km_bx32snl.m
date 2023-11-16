function bx32snl = km_bx32snl(s,n,l,csn,bsl,wps,lambdas,wcs,th)
% filename: km_bx32snl.m
% Calculate the coefficients of bx32snl for x-component of oblique 
% plasma wave model with a kappa-Maxweillian distrubution.
% Modified on Sep 24th, 2023

params_with_unit;

In = besseli(n,lambdas(s));

bx32snl = -1i*epsilon0*wps(s).^2*tan(th)*exp(-1*lambdas(s))/lambdas(s)...
          *n*In*bsl(s,l)*(csn(s,n)-n*wcs(s))/wcs(s);
end