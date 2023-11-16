function bx12snl = km_bx12snl(s,n,l,bsl,wps,lambdas)
% filename: km_bx12snl.m
% Calculate the coefficients of bx12snl for x-component of oblique 
% plasma waves with a kappa-Maxweillian distrubution.
% Modified on Sep 24th, 2023

params_with_unit;

In = besseli(n,lambdas(s));

bx12snl = -1i*epsilon0*wps(s).^2*exp(-1*lambdas(s))/lambdas(s)*n^2*In*bsl(s,l);
end