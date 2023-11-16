function M = M_bikappa(S,N,wcs,bxyzsn,MatrixNo_bikappa,ExyzNo)
% @Description: Compute the matrix M_bikappa for the perpendicular 
% propagation plasma wave model with a bi-kappa distribution.
% @Filename: M_bikappa.m
% @Date: 2021-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

SN = (2*N+1)*S;
index = getIndexOfBlkMatrix_bikappa(S,N,MatrixNo_bikappa);
FirstIndex = index(1)-1;

% create Matrix
len_row = SN+1;
len_col = 5*(SN+1) + 7;
M = zeros(len_row,len_col);

Nvector = -N:N;
for s=1:S
    for in=1:length(Nvector)
        n = Nvector(in);
        sn = (s-1)*length(Nvector) + in;
        M(sn,FirstIndex+sn) = M(sn,FirstIndex+sn) + n*wcs(s);
        M(sn,end-ExyzNo) = M(sn,end-ExyzNo) + bxyzsn(s,n);
    end
end

for s=1:S
    for in=1:length(Nvector)
        n = Nvector(in);
        sn = (s-1)*length(Nvector) + in;
        M(SN+1,FirstIndex+sn) = M(SN+1,FirstIndex+sn) + n*wcs(s);
        M(SN+1,end-ExyzNo) = M(SN+1,end-ExyzNo) + bxyzsn(s,n);
    end
end

