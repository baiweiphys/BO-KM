function bx30 = maxwell_bx30(S,N,J,bsnj,csnj,wps,lambdas,theta)
% filename: maxwell_bx30.m
% Calculate the coefficients of bx30 for the oblique plasma wave model 
% with a bi-Maxweillian distrubution.
% The script created by Bai Wei (email:baiwei12@mail.ustc.edu.cn) 
% was created on Sep 3rd, 2023，
% and revised on Sep 4th, 2023.

params_with_unit;
In = @(s,n) besseli(n,lambdas(s));

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

bx30 = -1i*epsilon0*tan(theta)*(coef1 + coef2);





