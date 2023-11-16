function index = getIndexOfBlkMatrix_bikappa(S,N,MatrixNo_bikappa)
% @Description: To obtain the index of a subblock matrix within a matrix 
% for perpendicular propagation in bi-kappa plasmas.
% @Filename: getIndexOfBlkMatrix_bikappa.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2021-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

SN = S*(2*N+1);
SNp1 = SN + 1;

switch MatrixNo_bikappa
    case 1
        index = 1:SNp1;
    case 2
        index = SNp1+1:2*SNp1;
    case 3
        index = 2*SNp1+1:3*SNp1;
    case 4
        index = 3*SNp1+1:4*SNp1;
    case 5
        index = 4*SNp1+1:5*SNp1;
    otherwise 
        disp('MatrixNo must be a integrer for in range 1 to 5.');
end

