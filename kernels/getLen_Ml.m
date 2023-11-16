function len = getLen_Ml(S,N,kappas)
% filename: getLen_Ml.m
% To get the length of Matrix Ml.
% Modified on Aug 11st, 2023

lmin = zeros(S,1);
lmax = zeros(S,1);
LL = zeros(S,1);

for s = 1:S
    lmin(s) = 1;
    lmax(s) = kappas(s);
    LL(s) = sum(lmin(s):lmax(s)); 
end
SNLL = (2*N+1)*sum(LL);
len = SNLL + 1;
