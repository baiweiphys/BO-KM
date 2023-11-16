function index = getIndexOfBlkMatrix_maxwell(S,N,J,MatrixNo_maxwell)
% filename: getIndexOfBlkMatrix_maxwell.m
% To get the index of subblock Matrix in total Matrix.
% The script created by Bai Wei (email:baiwei12@mail.ustc.edu.cn) 
% was created on Sep 3rd, 2023，
% and revised on Sep 4th, 2023.

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


