function by21snl = km_by21snl(s,n,l,bsnl,wps,lambdas)
% filename: km_by12snl.m
% Calculate the coefficients of by21snl for y-component of oblique 
% plasma wave model with loss-cone kappa-Maxweillian distrubution.
% Modified on Sep 24th, 2023

params_with_unit;

% dIn = besseli(n+1,lambdas(s)) + n*besseli(n,lambdas(s))/lambdas(s);
dIn = 0.5*besseli(n+1,lambdas(s)) + 0.5*besseli(n-1,lambdas(s));
In = besseli(n,lambdas(s));

by21snl = -1i*epsilon0*wps(s).^2*exp(-1*lambdas(s))/lambdas(s)*bsnl(s,n,l)...
          *(n^2*In + 2*lambdas(s)^2*(In-dIn));
end