function bx21snl = km_bx21snl(s,n,l,bsnl,wps,lambdas)
% filename: km_bx21snl.m
% Calculate the coefficients of bx21snl for x-component of oblique 
% plasma waves with a kappa-Maxweillian distrubution.
% Modified on Sep 24th, 2023

params_with_unit;

% dIn = besseli(n+1,lambdas(s)) + n*besseli(n,lambdas(s))/lambdas(s);
dIn = 0.5*besseli(n+1,lambdas(s)) + 0.5*besseli(n-1,lambdas(s));
In = besseli(n,lambdas(s));

bx21snl = -1*epsilon0*wps(s).^2*exp(-1*lambdas(s))*n*(In-dIn)*bsnl(s,n,l);
end