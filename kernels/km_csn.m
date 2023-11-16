function csn = km_csn(s,n,kappas,vts_parallel,wcs,us0,kz)
% filename: km_bx_KM.m
% Calculate the coefficients of csn(s) for x-component of oblique 
% plasma waves with a kappa-Maxweillian distrubution.
% Modified on Sep 24th, 2023

csn = n*wcs(s) + kz*us0(s) -1i*sqrt(kappas(s))*kz*vts_parallel(s);


