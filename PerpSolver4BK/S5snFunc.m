function S5sn = S5nFunc(s,n,kappas,lambdas)
% filename: S5snFunc.m
% Calculate the integral of S5sn, which is the expression derived from 
% the Summers1994 Pop.
% Generated on October 1st, 2023 by Bai Wei (baiweiphys@gmail.com).
% Modified on Sep 24th, 2023

dJn = @(x) -0.5*besselj(n+1,x) + 0.5*besselj(n-1,x);

tmp = @(x) x.^3.*dJn(x).^2./((1+0.5*x.^2/lambdas(s)/kappas(s)).^(kappas(s)+1.5)); 
tmp_int = integral(tmp,0,Inf);

S5sn = (1.0-0.25/kappas(s)^2)*tmp_int;
