function bsl = km_bsl(s,l,kappas,vts_parallel,vts_perp,kz)
% filename: km_bsl.m
% Calculate the coefficients for oblique plasma wave model with 
% a kappa-Maxweillian distrubution.
% Modified on Sep 24th, 2023

coefbs0 =  @(s,l) 1i/vts_parallel(s)*(kappas(s)-0.5)/sqrt(kappas(s));
coefbs1 =  @(s,l) coefbs0(s,l)*factorial(kappas(s))/factorial(2*kappas(s))...
          *factorial(2*kappas(s)-l-1)/factorial(kappas(s)-l);
coefbs2 = @(s,l) coefbs1(s,l)*(2i*sqrt(kappas(s))*kz*vts_parallel(s))^l;

bsl = coefbs2(s,l)*l*kz*vts_perp(s)^2;