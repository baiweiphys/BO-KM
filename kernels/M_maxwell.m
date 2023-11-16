function M = M_maxwell(S,N,J,csnj,bxyzsnj,MatrixNo_maxwell,ExyzNo)
% @Description: Calculate the M_maxwell matrix for the oblique 
% plasma wave model with a BM plasmas.
% @Filename: M_maxwell.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

NJ = (2*N+1)*J;
SNJ = S*NJ;

index = getIndexOfBlkMatrix_maxwell(S,N,J,MatrixNo_maxwell);
FirstIndex = index(1)-1;

% create Matrix
len_row = SNJ+1;
len_col = 9*(SNJ+1) + 9;
M = zeros(len_row,len_col);

Nvector = -N:N;
for s=1:S
    for in=1:length(Nvector)
        n = Nvector(in);
        for jj=1:J
            snj = (s-1)*NJ + (in-1)*J + jj;
            M(snj,FirstIndex+snj) = M(snj,FirstIndex+snj) + csnj(s,n,jj);
            M(snj,end-ExyzNo) = M(snj,end-ExyzNo) + bxyzsnj(s,n,jj);
        end
    end
end

for s=1:S
    for in=1:length(Nvector)
        n = Nvector(in);
        for jj=1:J
            snj = (s-1)*NJ + (in-1)*J + jj;
            M(SNJ+1,FirstIndex+snj) = M(SNJ+1,FirstIndex+snj) + csnj(s,n,jj);
            M(SNJ+1,end-ExyzNo) = M(SNJ+1,end-ExyzNo) + bxyzsnj(s,n,jj);
        end
    end
end


