function len = getLen_Mlp1(S,N,kappas)
% @Description: To get the size of Matrix Mlp1.
% @Filename: getLen_Mlp1.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-08-11
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

lmin = zeros(S,1);
lmax = zeros(S,1);
LL = zeros(S,1);

for s = 1:S
    lmin(s) = 1;
    lmax(s) = kappas(s);
    LL(s) = sum(lmin(s)+1:lmax(s)+1); 
end
SNLL = (2*N+1)*sum(LL); 
len = SNLL + 1;
