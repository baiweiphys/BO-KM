function by20 = maxwell_by20(S,N,J,bsnj,csnj,wps,lambdas)
% filename: maxwell_by20.m
% Calculate the coefficients of by20 for the oblique plasma wave model 
% with a bi-Maxweillian distrubution.
% The script created by Bai Wei (email:baiwei12@mail.ustc.edu.cn) 
% was created on Sep 3rd, 2023，
% and revised on Sep 4th, 2023.

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
                    *(n^2*In(s,n) + 2*lambdas(s)^2*(In(s,n)-dIn(s,n)))...
                    *bsnj(s,n,jj)/csnj(s,n,jj);
        end
    end
end

by20 = 1i*epsilon0*(coef1 + coef2);

