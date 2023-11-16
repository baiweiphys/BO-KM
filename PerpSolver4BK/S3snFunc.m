function S3sn = S3snFunc(s,n,kappas,lambdas)
% @Description: Compute the integral of S3sn, which is the expression 
% obtained from Ref. (Summers1994, PoP) .
% @Filename: S3snFunc.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

Jn = @(x) besselj(n,x);
dJn = @(x) -0.5*besselj(n+1,x) + 0.5*besselj(n-1,x);

tmp = @(x) x.^2.*Jn(x).*dJn(x)./((1+0.5*x.^2/lambdas(s)/kappas(s)).^(kappas(s)+1.5)); 
tmp_int = integral(tmp,0,Inf);

S3sn = (1.0-0.25/kappas(s)^2)*tmp_int;
