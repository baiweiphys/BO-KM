function [M,sum_b,index] = Mlm2_km(S,N,kappas,csn,bxyzsnl,MatrixNo,ExyzNo)
% @Description: To obtain the Mlm2_km matrix for KM plasmas.
% @Filename: Mlm2_km.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-08-12
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

% ExyzNo = 5 for Ex
% ExyzNo = 4 for Ey
% ExyzNo = 3 for Ez

sum_b = 0;

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
    LL(s) = sum(lmin(s):lmax(s)-2); 
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
        sum_b = sum_b - bxyzsnl(s,n,1);
        for l = lmin(s):lmax(s)
            if (l>2)
                for jj = 1:l-2
                    snll = (in-1)*sum(LL) + LLs(s) + sum(lmin(s):(l-3)) + jj;
                    M(snll,FirstIndex+snll) = M(snll,FirstIndex+snll) + csn(s,n);
                end
                if (jj<(l-2))  
                    M(snll,FirstIndex+snll+1) = M(snll,FirstIndex+snll+1) + 1;
                else
                    M(snll,end-ExyzNo) = M(snll,end-ExyzNo) + bxyzsnl(s,n,l);
                end
            end
        end
    end
end


% Step 4: jxyz
for in=1:length(Nvector)
    n = Nvector(in);
    for s=1:S
        for l = lmin(s):lmax(s) 
            if (l<2)
                M(SNLL+1,end-ExyzNo) = M(SNLL+1,end-ExyzNo) - csn(s,n)*bxyzsnl(s,n,1);
            elseif (l<3)
                M(SNLL+1,end-ExyzNo) = M(SNLL+1,end-ExyzNo) + bxyzsnl(s,n,l);
            else
                jj = 1; 
                nll = (in-1)*sum(LL) + LLs(s) + sum(lmin(s):(l-3)) + jj;
                M(SNLL+1,FirstIndex+snll) = M(SNLL+1,FirstIndex+snll) + 1;
            end
        end
    end
end
