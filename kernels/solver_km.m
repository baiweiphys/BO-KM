function w = solver_km(N,kx,kz,theta,B0,par,sp)
% @Description: To compute the roots for the oblique plasma 
% wave model with a kappa-Maxwellian distribution, given kx and kz.
% @Filename: solver_km.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-9-24
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

params_with_unit;

[S,index_skm,index_sbm,kappas,vts_parallel,vts_perp,Ts_parallel,Ts_perp,wps,wcs,us0,rhocs,lambdaDs] = getPlasmaPrameters(B0,par); 
lambdas = 0.5*kx^2*rhocs.^2; % for argument of the modified bessel function

csn = @(s,n) km_csn(s,n,kappas,vts_parallel,wcs,us0,kz);
bsl = @(s,l) km_bsl(s,l,kappas,vts_parallel,vts_perp,kz);
bsnl = @(s,n,l) km_bsnl(s,n,l,kappas,vts_parallel,wcs,kz);

% for x-component
bx10 = km_bx10(wps);
bx30 = km_bx30(wps,theta);
bx11snl = @(s,n,l) km_bx11snl(s,n,l,bsnl,wps,lambdas);
bx12snl = @(s,n,l) km_bx12snl(s,n,l,bsl,wps,lambdas);
bx21snl = @(s,n,l) km_bx21snl(s,n,l,bsnl,wps,lambdas);
bx22snl = @(s,n,l) km_bx22snl(s,n,l,bsl,wps,lambdas);
bx31snl = @(s,n,l) km_bx31snl(s,n,l,csn,bsl,bsnl,wps,lambdas,wcs,theta);
bx32snl = @(s,n,l) km_bx32snl(s,n,l,csn,bsl,wps,lambdas,wcs,theta);
bx33snl = @(s,n,l) km_bx33snl(s,n,l,bsnl,wps,lambdas,wcs,theta);

% for y-component
by11snl = @(s,n,l) -1*bx21snl(s,n,l);
by12snl = @(s,n,l) -1*bx22snl(s,n,l);
by20 = bx10;
by21snl = @(s,n,l) km_by21snl(s,n,l,bsnl,wps,lambdas);
by22snl = @(s,n,l) km_by22snl(s,n,l,bsl,wps,lambdas);
by31snl = @(s,n,l) km_by31snl(s,n,l,csn,bsl,bsnl,wps,lambdas,wcs,theta);
by32snl = @(s,n,l) km_by32snl(s,n,l,csn,bsl,wps,lambdas,wcs,theta);
by33snl = @(s,n,l) km_by33snl(s,n,l,bsnl,wps,lambdas,wcs,theta);


% for z-component
bz10 = bx30;
bz11snl = @(s,n,l) bx31snl(s,n,l);
bz12snl = @(s,n,l) bx32snl(s,n,l);
bz13snl = @(s,n,l) bx33snl(s,n,l);
bz21snl = @(s,n,l) -1*by31snl(s,n,l);
bz22snl = @(s,n,l) -1*by32snl(s,n,l);
bz23snl = @(s,n,l) -1*by33snl(s,n,l);

bz30 = km_bz30(wps,theta);
bz31snl = @(s,n,l) km_bz31snl(s,n,l,csn,bsl,bsnl,wps,lambdas,wcs,theta);
bz32snl = @(s,n,l) km_bz32snl(s,n,l,csn,bsl,wps,lambdas,wcs,theta);
bz33snl = @(s,n,l) km_bz33snl(s,n,l,csn,bsl,bsnl,wps,lambdas,wcs,theta);
bz34snl = @(s,n,l) km_bz34snl(s,n,l,bsnl,wps,lambdas,wcs,theta);

% Assembly matrix
%%%%%%%%%%%%%%%%%%%
% jxyzNo = 8 for jx
% jxyzNo = 7 for jy
% jxyzNo = 6 for jz
% ExyzNo = 5 for Ex
% ExyzNo = 4 for Ey
% ExyzNo = 3 for Ez
% BxyzNo = 2 for Bx
% BxyzNo = 1 for By
% BxyzNo = 0 for Bz
%%%%%%%%%%%%%%%%%%%

% Step 1
% for x-component
Mx11 = Ml_km(S,N,kappas,csn,bx11snl,1,5); % No.1
Mx21 = Ml_km(S,N,kappas,csn,bx21snl,2,4); % No.2
Mx31 = Ml_km(S,N,kappas,csn,bx31snl,3,3); % No.3
% 
Mx12 = Mlp1_km(S,N,kappas,csn,bx12snl,4,5); % No.4
Mx22 = Mlp1_km(S,N,kappas,csn,bx22snl,5,4); % No.5
Mx32 = Mlp1_km(S,N,kappas,csn,bx32snl,6,3); % No.6
% 
Mx33 = Mlm1_km(S,N,kappas,csn,bx33snl,7,3); % No.7
 
 
% for y-component
My11 = Ml_km(S,N,kappas,csn,by11snl,8,5); % No.8
My21 = Ml_km(S,N,kappas,csn,by21snl,9,4); % No.9
My31 = Ml_km(S,N,kappas,csn,by31snl,10,3); % No.10
% 
My12 = Mlp1_km(S,N,kappas,csn,by12snl,11,5); % No.11
My22 = Mlp1_km(S,N,kappas,csn,by22snl,12,4); % No.12
My32 = Mlp1_km(S,N,kappas,csn,by32snl,13,3); % No.13
% 
My33 = Mlm1_km(S,N,kappas,csn,by33snl,14,3); % No.14
 
% for z-component
Mz11 = Ml_km(S,N,kappas,csn,bz11snl,15,5); % No.15
Mz21 = Ml_km(S,N,kappas,csn,bz21snl,16,4); % No.16
Mz31 = Ml_km(S,N,kappas,csn,bz31snl,17,3); % No.17
% 
Mz12 = Mlp1_km(S,N,kappas,csn,bz12snl,18,5); % No.18
Mz22 = Mlp1_km(S,N,kappas,csn,bz22snl,19,4); % No.19
Mz32 = Mlp1_km(S,N,kappas,csn,bz32snl,20,3); % No.20
% 
Mz13 = Mlm1_km(S,N,kappas,csn,bz13snl,21,5); % No.21
Mz23 = Mlm1_km(S,N,kappas,csn,bz23snl,22,4); % No.22
Mz33 = Mlm1_km(S,N,kappas,csn,bz33snl,23,3); % No.23
% 
[Mz34,sum_b34,index34] = Mlm2_km(S,N,kappas,csn,bz34snl,24,3); % No.24

% Step 2
MatrixColLen = size(Mz34,2);
O = zeros(9,MatrixColLen);

M = [Mx11;Mx21;Mx31;Mx12;Mx22;Mx32;Mx33;...
     My11;My21;My31;My12;My22;My32;My33;...
     Mz11;Mz21;Mz31;Mz12;Mz22;Mz32;Mz13;Mz23;Mz33;Mz34;O];

A = eye(size(M)); 
A(index34(end),end-3) = sum_b34; % ExyzNo = 3 for Ez

% for jx
M(end-8,end-5) = M(end-8,end-5) + bx10; % for the first term of Ex
M(end-8,end-3) = M(end-8,end-3) + bx30; % for the second term of Ez

% for jy
M(end-7,end-4) = M(end-7,end-4) + by20; % for Ey

% for jz
M(end-6,end-5) = M(end-6,end-5) + bz10; % for the first term of Ex
M(end-6,end-3) = M(end-6,end-3) + bz30; % for the second term of Ez


% for Ex
M(end-5,end-1) = M(end-5,end-1) + c2*kz; % for the first item of Ex
endIndex_jx11 = getIndexOfBlkMatrix_km(S,N,kappas,1);
endIndex_jx21 = getIndexOfBlkMatrix_km(S,N,kappas,2);
endIndex_jx31 = getIndexOfBlkMatrix_km(S,N,kappas,3);
endIndex_jx12 = getIndexOfBlkMatrix_km(S,N,kappas,4);
endIndex_jx22 = getIndexOfBlkMatrix_km(S,N,kappas,5);
endIndex_jx32 = getIndexOfBlkMatrix_km(S,N,kappas,6);
endIndex_jx33 = getIndexOfBlkMatrix_km(S,N,kappas,7);
M(end-5,endIndex_jx11(end)) = M(end-5,endIndex_jx11(end)) - 1i/epsilon0; % for the second item of Ex
M(end-5,endIndex_jx21(end)) = M(end-5,endIndex_jx21(end)) - 1i/epsilon0; % for the second item of Ex
M(end-5,endIndex_jx31(end)) = M(end-5,endIndex_jx31(end)) - 1i/epsilon0; % for the second item of Ex
M(end-5,endIndex_jx12(end)) = M(end-5,endIndex_jx12(end)) - 1i/epsilon0; % for the second item of Ex
M(end-5,endIndex_jx22(end)) = M(end-5,endIndex_jx22(end)) - 1i/epsilon0; % for the second item of Ex
M(end-5,endIndex_jx32(end)) = M(end-5,endIndex_jx32(end)) - 1i/epsilon0; % for the second item of Ex
M(end-5,endIndex_jx33(end)) = M(end-5,endIndex_jx33(end)) - 1i/epsilon0; % for the second item of Ex
M(end-5,end-8) = M(end-5,end-8) - 1i/epsilon0; % for the second item of Ex with jx

% for Ey
M(end-4,end-2) = M(end-4,end-2) - c2*kz; % for the first item of Ey
M(end-4,end) = M(end-4,end) + c2*kx; % for the second item of Ey
endIndex_jy11 = getIndexOfBlkMatrix_km(S,N,kappas,8);
endIndex_jy21 = getIndexOfBlkMatrix_km(S,N,kappas,9);
endIndex_jy31 = getIndexOfBlkMatrix_km(S,N,kappas,10);
endIndex_jy12 = getIndexOfBlkMatrix_km(S,N,kappas,11);
endIndex_jy22 = getIndexOfBlkMatrix_km(S,N,kappas,12);
endIndex_jy32 = getIndexOfBlkMatrix_km(S,N,kappas,13);
endIndex_jy33 = getIndexOfBlkMatrix_km(S,N,kappas,14);
M(end-4,endIndex_jy11(end)) = M(end-4,endIndex_jy11(end)) - 1i/epsilon0; % for the third item of Ey
M(end-4,endIndex_jy21(end)) = M(end-4,endIndex_jy21(end)) - 1i/epsilon0; % for the third item of Ey
M(end-4,endIndex_jy31(end)) = M(end-4,endIndex_jy31(end)) - 1i/epsilon0; % for the third item of Ey
M(end-4,endIndex_jy12(end)) = M(end-4,endIndex_jy12(end)) - 1i/epsilon0; % for the third item of Ey
M(end-4,endIndex_jy22(end)) = M(end-4,endIndex_jy22(end)) - 1i/epsilon0; % for the third item of Ey
M(end-4,endIndex_jy32(end)) = M(end-4,endIndex_jy32(end)) - 1i/epsilon0; % for the third item of Ey
M(end-4,endIndex_jy33(end)) = M(end-4,endIndex_jy33(end)) - 1i/epsilon0; % for the third item of Ey
M(end-4,end-7) = M(end-4,end-7) - 1i/epsilon0; % for the third item of Ey


% Ez
M(end-3,end-1) = M(end-3,end-1) - c2*kx;  % for the first item of Ez
endIndex_jz11 = getIndexOfBlkMatrix_km(S,N,kappas,15);
endIndex_jz21 = getIndexOfBlkMatrix_km(S,N,kappas,16);
endIndex_jz31 = getIndexOfBlkMatrix_km(S,N,kappas,17);
endIndex_jz12 = getIndexOfBlkMatrix_km(S,N,kappas,18);
endIndex_jz22 = getIndexOfBlkMatrix_km(S,N,kappas,19);
endIndex_jz32 = getIndexOfBlkMatrix_km(S,N,kappas,20);
endIndex_jz13 = getIndexOfBlkMatrix_km(S,N,kappas,21);
endIndex_jz23 = getIndexOfBlkMatrix_km(S,N,kappas,22);
endIndex_jz33 = getIndexOfBlkMatrix_km(S,N,kappas,23);
endIndex_jz34 = getIndexOfBlkMatrix_km(S,N,kappas,24);
M(end-3,endIndex_jz11(end)) = M(end-3,endIndex_jz11(end)) - 1i/epsilon0; % for the second item of Ez
M(end-3,endIndex_jz21(end)) = M(end-3,endIndex_jz21(end)) - 1i/epsilon0; % for the second item of Ez
M(end-3,endIndex_jz31(end)) = M(end-3,endIndex_jz31(end)) - 1i/epsilon0; % for the second item of Ez
M(end-3,endIndex_jz12(end)) = M(end-3,endIndex_jz12(end)) - 1i/epsilon0; % for the second item of Ez
M(end-3,endIndex_jz22(end)) = M(end-3,endIndex_jz22(end)) - 1i/epsilon0; % for the second item of Ez
M(end-3,endIndex_jz32(end)) = M(end-3,endIndex_jz32(end)) - 1i/epsilon0; % for the second item of Ez
M(end-3,endIndex_jz13(end)) = M(end-3,endIndex_jz13(end)) - 1i/epsilon0; % for the second item of Ez
M(end-3,endIndex_jz23(end)) = M(end-3,endIndex_jz23(end)) - 1i/epsilon0; % for the second item of Ez
M(end-3,endIndex_jz33(end)) = M(end-3,endIndex_jz33(end)) - 1i/epsilon0; % for the second item of Ez
M(end-3,endIndex_jz34(end)) = M(end-3,endIndex_jz34(end)) - 1i/epsilon0; % for the second item of Ez
M(end-3,end-6) = M(end-3,end-6) - 1i/epsilon0; % for the second item of Ez with jz

% for Bx, By and Bz
M(end-2,end-4) = M(end-2,end-4) - kz;  % for Bx
M(end-1,end-5) = M(end-1,end-5) + kz; % for the first item of By
M(end-1,end-3) = M(end-1,end-3) - kx;      % for the second item of By
M(end,end-4) = M(end,end-4) + kx;   % for Bz


%% Solver
if(sp==0)
    [V,D] = eig(A\M);
else
    [V,D] = eigs(sparse(A\M),size(M,1),eps); 
end

omega = diag(D);
[wi,ind]=sort(imag(omega),'descend');
% [wr,ind]=sort(real(omega),'descend');
w = omega(ind);
w(abs(w)==min(abs(w))) = NaN+1i*NaN; % remove zero solution



