function S0sn = S0snFunc(s,n,kappas,lambdas)
% filename: S0snFunc.m
% Calculate the integral of S0sn, which is the expression derived from 
% the Summers1994 Pop.
% Generated on October 1st, 2023 by Bai Wei (baiweiphys@gmail.com).
% Modified on Sep 24th, 2023

Jn = @(x) besselj(n,x);

tmp = @(x) x.*Jn(x).^2./((1+0.5*x.^2/lambdas(s)/kappas(s)).^(kappas(s)+0.5)); 
tmp_int = integral(tmp,0,Inf);

S0sn = (1.0-0.5/kappas(s))*tmp_int;
