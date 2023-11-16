function index = getIndexOfBlkMatrix_maxwell(S,N,J,MatrixNo_maxwell)
% @Description: To obtain the index of a subblock matrix within a composite 
% matrix that exhibits a BM plasma distribution.
% @Filename: getIndexOfBlkMatrix_maxwell.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

SNJ = S*(2*N+1)*J;
SNJp1 = SNJ + 1;

switch MatrixNo_maxwell
    case 1
        index = 1:SNJp1;
    case 2
        index = SNJp1+1:2*SNJp1;
    case 3
        index = 2*SNJp1+1:3*SNJp1;
    case 4
        index = 3*SNJp1+1:4*SNJp1;
    case 5
        index = 4*SNJp1+1:5*SNJp1;
    case 6
        index = 5*SNJp1+1:6*SNJp1;
    case 7
        index = 6*SNJp1+1:7*SNJp1;
    case 8
        index = 7*SNJp1+1:8*SNJp1;
    case 9
        index = 8*SNJp1+1:9*SNJp1;
    otherwise 
        disp('MatrixNo must be a integrer for in range 1 to 9.');
end


