function bz31snl = km_bz31snl(s,n,l,csn,bsl,bsnl,wps,lambdas,wcs,th)
% filename: km_bz31snl.m
% Calculate the coefficients of bz31snl for z-component of oblique 
% plasma wave model with a kappa-Maxweillian distrubution.
% Modified on Sep 24th, 2023

params_with_unit;

% dIn = besseli(n+1,lambdas(s)) + n*besseli(n,lambdas(s))/lambdas(s);
In = besseli(n,lambdas(s));

bz31snl = -1i*epsilon0*wps(s).^2*tan(th)^2*exp(-1*lambdas(s))/lambdas(s)*In...
          *(bsnl(s,n,l)*(csn(s,n)-n*wcs(s))^2 + 2*bsl(s,l)*(csn(s,n)-n*wcs(s)))/wcs(s)^2;

end