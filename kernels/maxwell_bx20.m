function bx20 = maxwell_bx20(S,N,J,bsnj,csnj,wps,lambdas)
% @Description: Calculate the coefficients of bx20 for the oblique 
% plasma wave model with a bi-Maxwellian distribution.
% @Filename: maxwell_bx20.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

params_with_unit;
In = @(s,n) besseli(n,lambdas(s));
dIn = @(s,n) 0.5*besseli(n+1,lambdas(s)) + 0.5*besseli(n-1,lambdas(s));

coef = 0;
Nvector = -N:N;
for s=1:S
    for in=1:length(Nvector)
        n = Nvector(in);
        for jj=1:J
            coef = coef + wps(s)^2*exp(-1*lambdas(s))...
                    *n*(In(s,n)-dIn(s,n))*bsnj(s,n,jj)/csnj(s,n,jj);
        end
    end
end

bx20 = epsilon0*coef;





