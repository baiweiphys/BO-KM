function by22snl = km_by22snl(s,n,l,bsl,wps,lambdas)
% filename: coef_by22snl.m
% Calculate the coefficients of by22snl for y-component of oblique 
% plasma waves with a kappa-Maxweillian distrubution.
% Modified on Sep 24th, 2023

params_with_unit;

% dIn = besseli(n+1,lambdas(s)) + n*besseli(n,lambdas(s))/lambdas(s);
dIn = 0.5*besseli(n+1,lambdas(s)) + 0.5*besseli(n-1,lambdas(s));
In = besseli(n,lambdas(s));

by22snl = -1i*epsilon0*wps(s).^2*exp(-1*lambdas(s))/lambdas(s)*bsl(s,l)...
          *(n^2*In + 2*lambdas(s)^2*(In-dIn));
end