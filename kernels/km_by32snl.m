function by32snl = km_by32snl(s,n,l,csn,bsl,wps,lambdas,wcs,th)
% filename: km_by32snl.m
% Calculate the coefficients of by32snl for y-component of oblique 
% plasma wave model with a kappa-Maxweillian distrubution.
% Modified on Sep 24th, 2023

params_with_unit;

% dIn = besseli(n+1,lambdas(s)) + n*besseli(n,lambdas(s))/lambdas(s);
dIn = 0.5*besseli(n+1,lambdas(s)) + 0.5*besseli(n-1,lambdas(s));
In = besseli(n,lambdas(s));

by32snl = epsilon0*wps(s).^2*tan(th)*exp(-1*lambdas(s))*(In-dIn)...
          *bsl(s,l)*(csn(s,n)-n*wcs(s))/wcs(s);
end