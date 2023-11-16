function M = Mlm1_km(S,N,kappas,csn,bxyzsnl,MatrixNo,ExyzNo)
% filename: Mlm1_km.m
% To get the Matrix Ml_KM.
% Modified on Aug 12nd, 2023

% ExyzNo = 5 for Ex
% ExyzNo = 4 for Ey
% ExyzNo = 3 for Ez

% Step 0
len_Ml = getLen_Ml(S,N,kappas);
len_Mlp1 = getLen_Mlp1(S,N,kappas);
len_Mlm1 = getLen_Mlm1(S,N,kappas);
len_Mlm2 = getLen_Mlm2(S,N,kappas);

% Step 1
lmin = zeros(S,1);
lmax = zeros(S,1);
LL = zeros(S,1);
for s = 1:S
    lmin(s) = 1;
    lmax(s) = kappas(s);
    LL(s) = sum(lmin(s):lmax(s)-1); 
end
cumsum_LL = cumsum(LL);
LLs = zeros(S,1);
LLs(2:end) = cumsum_LL(1:end-1);
SNLL = (2*N+1)*sum(LL);


% Step 2: create Matrix
index = getIndexOfBlkMatrix_km(S,N,kappas,MatrixNo);
FirstIndex = index(1)-1;
len_row = SNLL+1;
len_col = 9*len_Ml + 9*len_Mlp1 + 5*len_Mlm1 + len_Mlm2 + 9;
M = zeros(len_row,len_col);

% Step 3
Nvector = -N:N;
for in=1:length(Nvector)
    n = Nvector(in);
    for s=1:S
        for l = lmin(s):lmax(s)
            if (l>1)
                for jj = 1:l-1
                    snll = (in-1)*sum(LL) + LLs(s) + sum(lmin(s):(l-2)) + jj;
                    M(snll,FirstIndex+snll) = M(snll,FirstIndex+snll) + csn(s,n); 
                    if (jj<l-1) 
                        M(snll,FirstIndex+snll+1) = M(snll,FirstIndex+snll+1) + 1;
                    else
                        M(snll,end-ExyzNo) = M(snll,end-ExyzNo) + bxyzsnl(s,n,l);
                    end
                end
            end
        end
    end
end


% step 4: jxyz 
for in=1:length(Nvector)
    n = Nvector(in);
    for s=1:S
        for l = lmin(s):lmax(s)
            if (l<2)
                M(SNLL+1,end-ExyzNo) = M(SNLL+1,end-ExyzNo) + bxyzsnl(s,n,1);
            else 
                jj = 1; 
                snll = (in-1)*sum(LL) + LLs(s) + sum(lmin(s):(l-2)) + jj;
                M(SNLL+1,FirstIndex+snll) = M(SNLL+1,FirstIndex+snll) + 1;
            end
        end
    end
end

