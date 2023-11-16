function len = getLen_Mlm2(S,N,kappas)
% @Description: To get the size of Matrix Mlm2.
% @Filename: getLen_Mlm2.m
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
    LL(s) = sum(lmin(s):lmax(s)-2); 
end
SNLL = (2*N+1)*sum(LL); 
len = SNLL + 1;