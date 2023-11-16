function bx22snl = km_bx22snl(s,n,l,bsl,wps,lambdas)
% filename: km_bx22snl.m
% Calculate the coefficients of bx22snl for x-component of oblique 
% plasma waves with a kappa-Maxweillian distrubution.
% Modified on Sep 24th, 2023

params_with_unit;

% dIn = besseli(n+1,lambdas(s)) + n*besseli(n,lambdas(s))/lambdas(s);
dIn = 0.5*besseli(n+1,lambdas(s)) + 0.5*besseli(n-1,lambdas(s));
In = besseli(n,lambdas(s));

bx22snl = -1*epsilon0*wps(s).^2*exp(-1*lambdas(s))*n*(In-dIn)*bsl(s,l);

end