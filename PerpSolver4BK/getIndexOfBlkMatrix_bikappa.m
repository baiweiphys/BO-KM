function index = getIndexOfBlkMatrix_bikappa(S,N,MatrixNo_bikappa)
% filename: getIndexOfBlkMatrix_bikappa.m
% To determine the index of a sub-block matrix within a total matrix 
% for perpendicular propagation in bi-kappa plasmas.
% Generated on October 1st, 2023 by Bai Wei (baiweiphys@gmail.com).
% Modified on Sep 24th, 2023

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

