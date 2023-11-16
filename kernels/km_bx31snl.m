function bx31snl = km_bx31snl(s,n,l,csn,bsl,bsnl,wps,lambdas,wcs,th)
% filename: km_bx31snl.m
% Calculate the coefficients of bx31snl for x-component of oblique 
% plasma wave model with a kappa-Maxweillian distrubution.
% Modified on Sep 24th, 2023

params_with_unit;

In = besseli(n,lambdas(s));

bx31snl = -1i*epsilon0*wps(s).^2*tan(th)*exp(-1*lambdas(s))/lambdas(s)...
          *n*In*(bsnl(s,n,l)*(csn(s,n)-n*wcs(s))+bsl(s,l))/wcs(s);
end