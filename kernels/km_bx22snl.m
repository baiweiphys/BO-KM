function bx22snl = km_bx22snl(s,n,l,bsl,wps,lambdas)
% @Description: Calculate the coefficients of bx22snl for the x-component 
% of oblique plasma waves with a kappa-Maxwellian distribution.
% @Filename: km_bx22snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

params_with_unit;

% dIn = besseli(n+1,lambdas(s)) + n*besseli(n,lambdas(s))/lambdas(s);
dIn = 0.5*besseli(n+1,lambdas(s)) + 0.5*besseli(n-1,lambdas(s));
In = besseli(n,lambdas(s));

bx22snl = -1*epsilon0*wps(s).^2*exp(-1*lambdas(s))*n*(In-dIn)*bsl(s,l);

end