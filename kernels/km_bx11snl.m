function bx11snl = km_bx11snl(s,n,l,bsnl,wps,lambdas)
% filename: km_bx11snl.m
% Calculate the coefficients of bx11snl for x-component of oblique 
% plasma waves with a kappa-Maxweillian distrubution.
% Modified on Sep 24th, 2023

params_with_unit;

In = besseli(n,lambdas(s));

bx11snl = -1i*epsilon0*wps(s).^2*exp(-1*lambdas(s))/lambdas(s)*n^2*In*bsnl(s,n,l);
end