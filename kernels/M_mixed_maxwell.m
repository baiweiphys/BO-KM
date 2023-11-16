function M = M_mixed_maxwell(S_km,S_bm,N,J,kappas_km,csnj_maxwell,bxyzsnj,MatrixNo,ExyzNo)
% @Description: Calculate the M_mixed_maxwell matrix for the oblique 
% plasma wave model incorporating a hybrid distribution of KM and BM plasmas.
% @Filename: M_mixed_maxwell.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15


% Step 0
len_Ml = getLen_Ml(S_km,N,kappas_km);
len_Mlp1 = getLen_Mlp1(S_km,N,kappas_km);
len_Mlm1 = getLen_Mlm1(S_km,N,kappas_km);
len_Mlm2 = getLen_Mlm2(S_km,N,kappas_km);


% Step 2: get the length of the maxwellian matrix
index = getIndexOfBlkMatrix_mixed(S_km,S_bm,N,J,kappas_km,MatrixNo);
FirstIndex = index(1)-1;

% for kappa-Maxwellian matrix 
len_km_col = 9*len_Ml + 9*len_Mlp1 + 5*len_Mlm1 + len_Mlm2;
% for bi-Maxwellian matrix 
SNJ = S_bm*(2*N+1)*J;
len_maxwell_col = 9*(SNJ + 1);

len_row = SNJ+1;
len_col = len_km_col + len_maxwell_col + 9;
M = zeros(len_row,len_col);


% Step3: create Matrix
NJ = (2*N+1)*J;
SNJ = S_bm*NJ;

Nvector = -N:N;
for s=1:S_bm
    for in=1:length(Nvector)
        n = Nvector(in);
        for jj=1:J
            snj = (s-1)*NJ + (in-1)*J + jj;
            M(snj,FirstIndex+snj) = M(snj,FirstIndex+snj) + csnj_maxwell(s,n,jj);
            M(snj,end-ExyzNo) = M(snj,end-ExyzNo) + bxyzsnj(s,n,jj);
        end
    end
end

for s=1:S_bm
    for in=1:length(Nvector)
        n = Nvector(in);
        for jj=1:J
            snj = (s-1)*NJ + (in-1)*J + jj;
            M(SNJ+1,FirstIndex+snj) = M(SNJ+1,FirstIndex+snj) + csnj_maxwell(s,n,jj);
            M(SNJ+1,end-ExyzNo) = M(SNJ+1,end-ExyzNo) + bxyzsnj(s,n,jj);
        end
    end
end


