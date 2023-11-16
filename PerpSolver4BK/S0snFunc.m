function S0sn = S0snFunc(s,n,kappas,lambdas)
% @Description: Compute the integral of S0sn, which is the expression 
% obtained from Ref. (Summers1994, PoP) .
% @Filename: S0snFunc.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

Jn = @(x) besselj(n,x);

tmp = @(x) x.*Jn(x).^2./((1+0.5*x.^2/lambdas(s)/kappas(s)).^(kappas(s)+0.5)); 
tmp_int = integral(tmp,0,Inf);

S0sn = (1.0-0.5/kappas(s))*tmp_int;
