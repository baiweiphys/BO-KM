function by33snl = km_by33snl(s,n,l,bsnl,wps,lambdas,wcs,th)
% filename: km_by32snl.m
% Calculate the coefficients of by33snl for y-component of oblique 
% plasma wave model with a kappa-Maxweillian distrubution.
% Modified on Sep 24th, 2023

params_with_unit;

% dIn = besseli(n+1,lambdas(s)) + n*besseli(n,lambdas(s))/lambdas(s);
dIn = 0.5*besseli(n+1,lambdas(s)) + 0.5*besseli(n-1,lambdas(s));
In = besseli(n,lambdas(s));

by33snl = epsilon0*wps(s).^2*tan(th)*exp(-1*lambdas(s))*(In-dIn)...
          *bsnl(s,n,l)/wcs(s);
end