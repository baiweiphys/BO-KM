function index = getIndexOfBlkMatrix_km(S,N,kappas,MatrixNo)
% filename: getIndexOfBlkMatrix.m
% To get the index of KM subblock Matrix in total Matrix.
% Modified on Aug 11st, 2023

len_Ml = getLen_Ml(S,N,kappas);
len_Mlp1 = getLen_Mlp1(S,N,kappas);
len_Mlm1 = getLen_Mlm1(S,N,kappas);
len_Mlm2 = getLen_Mlm2(S,N,kappas);

switch MatrixNo
    case 1
        index = 1:len_Ml;
    case 2
        index = len_Ml+1:2*len_Ml;
    case 3
        index = 2*len_Ml+1:3*len_Ml;
    case 4
        index = 3*len_Ml+1:3*len_Ml+len_Mlp1;
    case 5
        index = 3*len_Ml+len_Mlp1+1:3*len_Ml+2*len_Mlp1;
    case 6
        index = 3*len_Ml+2*len_Mlp1+1:3*len_Ml+3*len_Mlp1;
    case 7
        index = 3*len_Ml+3*len_Mlp1+1:3*len_Ml+3*len_Mlp1+len_Mlm1;
    case 8
        index = 3*len_Ml+3*len_Mlp1+len_Mlm1+1:4*len_Ml+3*len_Mlp1+len_Mlm1;
    case 9
        index = 4*len_Ml+3*len_Mlp1+len_Mlm1+1:5*len_Ml+3*len_Mlp1+len_Mlm1;
    case 10
        index = 5*len_Ml+3*len_Mlp1+len_Mlm1+1:6*len_Ml+3*len_Mlp1+len_Mlm1;
    case 11
        index =6*len_Ml+3*len_Mlp1+len_Mlm1+1:6*len_Ml+4*len_Mlp1+len_Mlm1;
    case 12
        index =6*len_Ml+4*len_Mlp1+len_Mlm1+1:6*len_Ml+5*len_Mlp1+len_Mlm1;
    case 13
        index =6*len_Ml+5*len_Mlp1+len_Mlm1+1:6*len_Ml+6*len_Mlp1+len_Mlm1;
    case 14
        index =6*len_Ml+6*len_Mlp1+len_Mlm1+1:6*len_Ml+6*len_Mlp1+2*len_Mlm1;
    case 15
        index =6*len_Ml+6*len_Mlp1+2*len_Mlm1+1:7*len_Ml+6*len_Mlp1+2*len_Mlm1;
    case 16
        index =7*len_Ml+6*len_Mlp1+2*len_Mlm1+1:8*len_Ml+6*len_Mlp1+2*len_Mlm1;
    case 17
        index =8*len_Ml+6*len_Mlp1+2*len_Mlm1+1:9*len_Ml+6*len_Mlp1+2*len_Mlm1;
    case 18
        index =9*len_Ml+6*len_Mlp1+2*len_Mlm1+1:9*len_Ml+7*len_Mlp1+2*len_Mlm1;
    case 19
        index =9*len_Ml+7*len_Mlp1+2*len_Mlm1+1:9*len_Ml+8*len_Mlp1+2*len_Mlm1;
    case 20
        index =9*len_Ml+8*len_Mlp1+2*len_Mlm1+1:9*len_Ml+9*len_Mlp1+2*len_Mlm1;
    case 21
        index =9*len_Ml+9*len_Mlp1+2*len_Mlm1+1:9*len_Ml+9*len_Mlp1+3*len_Mlm1;
    case 22
        index =9*len_Ml+9*len_Mlp1+3*len_Mlm1+1:9*len_Ml+9*len_Mlp1+4*len_Mlm1;
    case 23
        index =9*len_Ml+9*len_Mlp1+4*len_Mlm1+1:9*len_Ml+9*len_Mlp1+5*len_Mlm1;
    case 24
        index =9*len_Ml+9*len_Mlp1+5*len_Mlm1+1:9*len_Ml+9*len_Mlp1+5*len_Mlm1+len_Mlm2;
    otherwise 
        disp('MatrixNo must be a integrer for in range 1 to 24.');
end


