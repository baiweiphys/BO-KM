function by31snl =km_by31snl(s,n,l,csn,bsl,bsnl,wps,lambdas,wcs,th)
% filename: km_by31snl.m
% Calculate the coefficients of by31snl for y-component of oblique 
% plasma wave model with loss-cone kappa-Maxweillian distrubution.
% Modified on Sep 24th, 2023

params_with_unit;

% dIn = besseli(n+1,lambdas(s)) + n*besseli(n,lambdas(s))/lambdas(s);
dIn = 0.5*besseli(n+1,lambdas(s)) + 0.5*besseli(n-1,lambdas(s));
In = besseli(n,lambdas(s));

by31snl = epsilon0*wps(s).^2*tan(th)*exp(-1*lambdas(s))*(In-dIn)...
          *(bsnl(s,n,l)*(csn(s,n)-n*wcs(s)) + bsl(s,l))/wcs(s);
end