function len = getLen_Mlm2(S,N,kappas)
% filename: getLen_Mlm2.m
% To get the length of Matrix Mlm2.
% Modified on Aug 11st, 2023

lmin = zeros(S,1);
lmax = zeros(S,1);
LL = zeros(S,1);

for s = 1:S
    lmin(s) = 1;
    lmax(s) = kappas(s);
    LL(s) = sum(lmin(s):lmax(s)-2); 
end
SNLL = (2*N+1)*sum(LL); 
len = SNLL + 1;