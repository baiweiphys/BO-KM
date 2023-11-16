function bz30 = maxwell_bz30(S,N,J,bsnj,csnj,wps,lambdas,theta)
% @Description: Calculate the coefficients of bz30 for the oblique 
% plasma wave model with a bi-Maxwellian distribution.
% @Filename: maxwell_bz30.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

params_with_unit;
In = @(s,n) besseli(n,lambdas(s));
dIn = @(s,n) 0.5*besseli(n+1,lambdas(s)) + 0.5*besseli(n-1,lambdas(s));

coef1 = sum(wps.^2);
coef2 = 0;
Nvector = -N:N;
for s=1:S
    for in=1:length(Nvector)
        n = Nvector(in);
        for jj=1:J
            coef2 = coef2 + wps(s)^2*exp(-1*lambdas(s))/lambdas(s)...
                    *n^2*In(s,n)*bsnj(s,n,jj)/csnj(s,n,jj);
        end
    end
end

bz30 = 1i*epsilon0*tan(theta)^2*(coef1 + coef2);
