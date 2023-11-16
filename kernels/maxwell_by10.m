function by10 = maxwell_by10(S,N,J,bsnj,csnj,wps,lambdas)
% filename: maxwell_by10.m
% Calculate the coefficients of by10 for the oblique plasma wave model 
% with a bi-Maxweillian distrubution.
% The script created by Bai Wei (email:baiwei12@mail.ustc.edu.cn) 
% was created on Sep 3rd, 2023，
% and revised on Sep 4th, 2023.

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

by10 = -1*epsilon0*coef;





