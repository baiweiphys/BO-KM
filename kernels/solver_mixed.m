function [w,dEx,dEy,dEz,dBx,dBy,dBz] = solver_mixed(N,J,kx,kz,theta,B0,par,sp)
% filename: solver_mixed.m
% To calcuate the roots by given k for oblique plasma wave model 
% with a mixed distrubutions in both kappa-Maxweillian and bi-Maxweillian plasmas.
% Modified on Oct 16th, 2023

params_with_unit;

[S,index_skm,index_sbm,kappas,vts_parallel,vts_perp,Ts_parallel,Ts_perp,wps,wcs,us0,rhocs,lambdaDs] = getPlasmaPrameters(B0,par);
lambdas = 0.5*kx^2*rhocs.^2; % for argument of the modified bessel function

% for kappa-Maxwellian 
S_km = sum(index_skm);
kappas_km = kappas(index_skm);
vts_parallel_km = vts_parallel(index_skm);
vts_perp_km = vts_perp(index_skm);
Ts_parallel_km = Ts_parallel(index_skm);
Ts_perp_km = Ts_perp(index_skm);
wps_km = wps(index_skm);
wcs_km = wcs(index_skm);
us0_km = us0(index_skm);
rhocs_km = rhocs(index_skm);
lambdaDs_km = lambdaDs(index_skm);
lambdas_km = lambdas(index_skm);

% for bi-Maxwellian plasma
S_bm = sum(index_sbm);
vts_parallel_bm = vts_parallel(index_sbm);
vts_perp_bm = vts_perp(index_sbm);
Ts_parallel_bm = Ts_parallel(index_sbm);
Ts_perp_bm = Ts_perp(index_sbm);
wps_bm = wps(index_sbm);
wcs_bm = wcs(index_sbm);
us0_bm = us0(index_sbm);
rhocs_bm = rhocs(index_sbm);
lambdaDs_bm = lambdaDs(index_sbm);
lambdas_bm = lambdas(index_sbm);

%% for kappa-Maxwellian plasma
csn_km = @(s,n) km_csn(s,n,kappas_km,vts_parallel_km,wcs_km,us0_km,kz);
bsl_km = @(s,l) km_bsl(s,l,kappas_km,vts_parallel_km,vts_perp_km,kz);
bsnl_km = @(s,n,l) km_bsnl(s,n,l,kappas_km,vts_parallel_km,wcs_km,kz);

% for x-component
bx10_km = km_bx10(wps_km);
bx30_km = km_bx30(wps_km,theta);
bx11snl = @(s,n,l) km_bx11snl(s,n,l,bsnl_km,wps_km,lambdas_km);
bx12snl = @(s,n,l) km_bx12snl(s,n,l,bsl_km,wps_km,lambdas_km);
bx21snl = @(s,n,l) km_bx21snl(s,n,l,bsnl_km,wps_km,lambdas_km);
bx22snl = @(s,n,l) km_bx22snl(s,n,l,bsl_km,wps_km,lambdas_km);
bx31snl = @(s,n,l) km_bx31snl(s,n,l,csn_km,bsl_km,bsnl_km,wps_km,lambdas_km,wcs_km,theta);
bx32snl = @(s,n,l) km_bx32snl(s,n,l,csn_km,bsl_km,wps_km,lambdas_km,wcs_km,theta);
bx33snl = @(s,n,l) km_bx33snl(s,n,l,bsnl_km,wps_km,lambdas_km,wcs_km,theta);

% for y-component
by11snl = @(s,n,l) -1*bx21snl(s,n,l);
by12snl = @(s,n,l) -1*bx22snl(s,n,l);
by20_km = bx10_km;
by21snl = @(s,n,l) km_by21snl(s,n,l,bsnl_km,wps_km,lambdas_km);
by22snl = @(s,n,l) km_by22snl(s,n,l,bsl_km,wps_km,lambdas_km);
by31snl = @(s,n,l) km_by31snl(s,n,l,csn_km,bsl_km,bsnl_km,wps_km,lambdas_km,wcs_km,theta);
by32snl = @(s,n,l) km_by32snl(s,n,l,csn_km,bsl_km,wps_km,lambdas_km,wcs_km,theta);
by33snl = @(s,n,l) km_by33snl(s,n,l,bsnl_km,wps_km,lambdas_km,wcs_km,theta);

% for z-component
bz10_km = bx30_km;
bz11snl = @(s,n,l) bx31snl(s,n,l);
bz12snl = @(s,n,l) bx32snl(s,n,l);
bz13snl = @(s,n,l) bx33snl(s,n,l);
bz21snl = @(s,n,l) -1*by31snl(s,n,l);
bz22snl = @(s,n,l) -1*by32snl(s,n,l);
bz23snl = @(s,n,l) -1*by33snl(s,n,l);

bz30_km = km_bz30(wps_km,theta);
bz31snl = @(s,n,l) km_bz31snl(s,n,l,csn_km,bsl_km,bsnl_km,wps_km,lambdas_km,wcs_km,theta);
bz32snl = @(s,n,l) km_bz32snl(s,n,l,csn_km,bsl_km,wps_km,lambdas_km,wcs_km,theta);
bz33snl = @(s,n,l) km_bz33snl(s,n,l,csn_km,bsl_km,bsnl_km,wps_km,lambdas_km,wcs_km,theta);
bz34snl = @(s,n,l) km_bz34snl(s,n,l,bsnl_km,wps_km,lambdas_km,wcs_km,theta);

%% for bi-Maxwellian plasma
bsnj_maxwell = @(s,n,jj) maxwell_bsnj(s,n,jj,vts_parallel_bm,Ts_parallel_bm,Ts_perp_bm,wcs_bm,kz,J);
csnj_maxwell = @(s,n,jj) maxwell_csnj(s,n,jj,vts_parallel_bm,wcs_bm,us0_bm,kz,J);

% for x-component
bx10_maxwell = maxwell_bx10(S_bm,N,J,bsnj_maxwell,csnj_maxwell,wps_bm,lambdas_bm);
bx11snj = @(s,n,jj) maxwell_bx11snj(s,n,jj,bsnj_maxwell,csnj_maxwell,wps_bm,lambdas_bm);
bx20_maxwell = maxwell_bx20(S_bm,N,J,bsnj_maxwell,csnj_maxwell,wps_bm,lambdas_bm);
bx21snj = @(s,n,jj) maxwell_bx21snj(s,n,jj,bsnj_maxwell,csnj_maxwell,wps_bm,lambdas_bm);
bx30_maxwell = maxwell_bx30(S_bm,N,J,bsnj_maxwell,csnj_maxwell,wps_bm,lambdas_bm,theta);
bx31snj = @(s,n,jj) maxwell_bx31snj(s,n,jj,bsnj_maxwell,csnj_maxwell,wps_bm,wcs_bm,lambdas_bm,theta);

% for y-component
by10_maxwell = maxwell_by10(S_bm,N,J,bsnj_maxwell,csnj_maxwell,wps_bm,lambdas_bm);
by11snj = @(s,n,jj) maxwell_by11snj(s,n,jj,bsnj_maxwell,csnj_maxwell,wps_bm,lambdas_bm);
by20_maxwell = maxwell_by20(S_bm,N,J,bsnj_maxwell,csnj_maxwell,wps_bm,lambdas_bm);
by21snj = @(s,n,jj) maxwell_by21snj(s,n,jj,bsnj_maxwell,csnj_maxwell,wps_bm,lambdas_bm);
by30_maxwell = maxwell_by30(S_bm,N,J,bsnj_maxwell,csnj_maxwell,wps_bm,lambdas_bm,theta);
by31snj = @(s,n,jj) maxwell_by31snj(s,n,jj,bsnj_maxwell,csnj_maxwell,wps_bm,wcs_bm,lambdas_bm,theta);

% for z-component
bz10_maxwell = maxwell_bz10(S_bm,N,J,bsnj_maxwell,csnj_maxwell,wps_bm,lambdas_bm,theta);
bz11snj = @(s,n,jj) maxwell_bz11snj(s,n,jj,bsnj_maxwell,csnj_maxwell,wps_bm,wcs_bm,lambdas_bm,theta);
bz20_maxwell = maxwell_bz20(S_bm,N,J,bsnj_maxwell,csnj_maxwell,wps_bm,lambdas_bm,theta);
bz21snj = @(s,n,jj) maxwell_bz21snj(s,n,jj,bsnj_maxwell,csnj_maxwell,wps_bm,wcs_bm,lambdas_bm,theta);
bz30_maxwell = maxwell_bz30(S_bm,N,J,bsnj_maxwell,csnj_maxwell,wps_bm,lambdas_bm,theta);
bz31snj = @(s,n,jj) maxwell_bz31snj(s,n,jj,bsnj_maxwell,csnj_maxwell,wps_bm,wcs_bm,lambdas_bm,theta);

%% Assemble matrix
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

%% Step 1: kappa-Maxwellian plasma matrix
% for x-component
Mx11_km = Ml_mixed_km(S_km,S_bm,N,J,kappas_km,csn_km,bx11snl,1,5); % No.1
Mx21_km = Ml_mixed_km(S_km,S_bm,N,J,kappas_km,csn_km,bx21snl,2,4); % No.2
Mx31_km = Ml_mixed_km(S_km,S_bm,N,J,kappas_km,csn_km,bx31snl,3,3); % No.3
%
Mx12_km = Mlp1_mixed_km(S_km,S_bm,N,J,kappas_km,csn_km,bx12snl,4,5); % No.4
Mx22_km = Mlp1_mixed_km(S_km,S_bm,N,J,kappas_km,csn_km,bx22snl,5,4); % No.5
Mx32_km = Mlp1_mixed_km(S_km,S_bm,N,J,kappas_km,csn_km,bx32snl,6,3); % No.6
% 
Mx33_km = Mlm1_mixed_km(S_km,S_bm,N,J,kappas_km,csn_km,bx33snl,7,3); % No.7 

% for y-component
My11_km = Ml_mixed_km(S_km,S_bm,N,J,kappas_km,csn_km,by11snl,8,5); % No.8
My21_km = Ml_mixed_km(S_km,S_bm,N,J,kappas_km,csn_km,by21snl,9,4); % No.9
My31_km = Ml_mixed_km(S_km,S_bm,N,J,kappas_km,csn_km,by31snl,10,3); % No.10
% 
My12_km = Mlp1_mixed_km(S_km,S_bm,N,J,kappas_km,csn_km,by12snl,11,5); % No.11
My22_km = Mlp1_mixed_km(S_km,S_bm,N,J,kappas_km,csn_km,by22snl,12,4); % No.12
My32_km = Mlp1_mixed_km(S_km,S_bm,N,J,kappas_km,csn_km,by32snl,13,3); % No.13
% 
My33_km = Mlm1_mixed_km(S_km,S_bm,N,J,kappas_km,csn_km,by33snl,14,3); % No.14 
 
% for z-component
Mz11_km = Ml_mixed_km(S_km,S_bm,N,J,kappas_km,csn_km,bz11snl,15,5); % No.15
Mz21_km = Ml_mixed_km(S_km,S_bm,N,J,kappas_km,csn_km,bz21snl,16,4); % No.16
Mz31_km = Ml_mixed_km(S_km,S_bm,N,J,kappas_km,csn_km,bz31snl,17,3); % No.17
% 
Mz12_km = Mlp1_mixed_km(S_km,S_bm,N,J,kappas_km,csn_km,bz12snl,18,5); % No.18
Mz22_km = Mlp1_mixed_km(S_km,S_bm,N,J,kappas_km,csn_km,bz22snl,19,4); % No.19
Mz32_km = Mlp1_mixed_km(S_km,S_bm,N,J,kappas_km,csn_km,bz32snl,20,3); % No.20
% 
Mz13_km = Mlm1_mixed_km(S_km,S_bm,N,J,kappas_km,csn_km,bz13snl,21,5); % No.21 
Mz23_km = Mlm1_mixed_km(S_km,S_bm,N,J,kappas_km,csn_km,bz23snl,22,4); % No.22 
Mz33_km = Mlm1_mixed_km(S_km,S_bm,N,J,kappas_km,csn_km,bz33snl,23,3); % No.23 
% 
[Mz34_km,sum_b34_km,index34_km] = Mlm2_mixed_km(S_km,S_bm,N,J,kappas_km,csn_km,bz34snl,24,3); % No.24


%% Step 2: bi-Maxwellian plasma matrix
% for x-component
Mx11_maxwell = M_mixed_maxwell(S_km,S_bm,N,J,kappas_km,csnj_maxwell,bx11snj,25,5); % No.25
Mx21_maxwell = M_mixed_maxwell(S_km,S_bm,N,J,kappas_km,csnj_maxwell,bx21snj,26,4); % No.26
Mx31_maxwell = M_mixed_maxwell(S_km,S_bm,N,J,kappas_km,csnj_maxwell,bx31snj,27,3); % No.27

% for y-component
My11_maxwell = M_mixed_maxwell(S_km,S_bm,N,J,kappas_km,csnj_maxwell,by11snj,28,5); % No.28
My21_maxwell = M_mixed_maxwell(S_km,S_bm,N,J,kappas_km,csnj_maxwell,by21snj,29,4); % No.29
My31_maxwell = M_mixed_maxwell(S_km,S_bm,N,J,kappas_km,csnj_maxwell,by31snj,30,3); % No.30

% for z-component
Mz11_maxwell = M_mixed_maxwell(S_km,S_bm,N,J,kappas_km,csnj_maxwell,bz11snj,31,5); % No.31
Mz21_maxwell = M_mixed_maxwell(S_km,S_bm,N,J,kappas_km,csnj_maxwell,bz21snj,32,4); % No.32
Mz31_maxwell = M_mixed_maxwell(S_km,S_bm,N,J,kappas_km,csnj_maxwell,bz31snj,33,3); % No.33

%% Step 3: Assemble matrix based on the mixed distributions of KM and BM plasmas.

M_km = [Mx11_km;Mx21_km;Mx31_km;Mx12_km;Mx22_km;Mx32_km;Mx33_km;...
        My11_km;My21_km;My31_km;My12_km;My22_km;My32_km;My33_km;...
        Mz11_km;Mz21_km;Mz31_km;Mz12_km;Mz22_km;Mz32_km;Mz13_km;Mz23_km;Mz33_km;Mz34_km];

M_maxwell = [Mx11_maxwell;Mx21_maxwell;Mx31_maxwell;...
             My11_maxwell;My21_maxwell;My31_maxwell;...
             Mz11_maxwell;Mz21_maxwell;Mz31_maxwell];

MatrixColLen_mixed = size(Mx11_maxwell,2);
O = zeros(9,MatrixColLen_mixed);

M = [M_km;M_maxwell;O];

A = eye(size(M)); 
A(index34_km(end),end-3) = sum_b34_km; % ExyzNo = 3 for Ez

%% Step 4: For the last equations of No.34, No.35, and No.36.
% for jx (No.34)
M(end-8,end-5) = M(end-8,end-5) + bx10_km + bx10_maxwell; % for the first term of Ex
M(end-8,end-3) = M(end-8,end-3) + bx30_km + bx30_maxwell; % for the second term of Ez
M(end-8,end-4) = M(end-8,end-4) + bx20_maxwell; % for the second term of Ey

% for jy (No.35)
M(end-7,end-5) = M(end-7,end-5) + by10_maxwell; % for the first term of Ex
M(end-7,end-4) = M(end-7,end-4) + by20_km + by20_maxwell; % for the second term of Ey
M(end-7,end-3) = M(end-7,end-3) + by30_maxwell; % for the third term of Ez

% for jz (No.36)
M(end-6,end-5) = M(end-6,end-5) + bz10_km + bz10_maxwell; % for the first term of Ex
M(end-6,end-4) = M(end-6,end-4) + bz20_maxwell; % for the second term of Ey
M(end-6,end-3) = M(end-6,end-3) + bz30_km + bz30_maxwell; % for the second term of Ez

%% Step 5: The perturbed currents of Maxwell's equations for the kappa-Maxwellian plasmas
% for Ex                   
Index_jx11_km = getIndexOfBlkMatrix_mixed(S_km,S_bm,N,J,kappas_km,1);
Index_jx21_km = getIndexOfBlkMatrix_mixed(S_km,S_bm,N,J,kappas_km,2);
Index_jx31_km = getIndexOfBlkMatrix_mixed(S_km,S_bm,N,J,kappas_km,3);
Index_jx12_km = getIndexOfBlkMatrix_mixed(S_km,S_bm,N,J,kappas_km,4);
Index_jx22_km = getIndexOfBlkMatrix_mixed(S_km,S_bm,N,J,kappas_km,5);
Index_jx32_km = getIndexOfBlkMatrix_mixed(S_km,S_bm,N,J,kappas_km,6);
Index_jx33_km = getIndexOfBlkMatrix_mixed(S_km,S_bm,N,J,kappas_km,7);
M(end-5,Index_jx11_km(end)) = M(end-5,Index_jx11_km(end)) - 1i/epsilon0; % for the second item of Ex
M(end-5,Index_jx21_km(end)) = M(end-5,Index_jx21_km(end)) - 1i/epsilon0; % for the second item of Ex
M(end-5,Index_jx31_km(end)) = M(end-5,Index_jx31_km(end)) - 1i/epsilon0; % for the second item of Ex
M(end-5,Index_jx12_km(end)) = M(end-5,Index_jx12_km(end)) - 1i/epsilon0; % for the second item of Ex
M(end-5,Index_jx22_km(end)) = M(end-5,Index_jx22_km(end)) - 1i/epsilon0; % for the second item of Ex
M(end-5,Index_jx32_km(end)) = M(end-5,Index_jx32_km(end)) - 1i/epsilon0; % for the second item of Ex
M(end-5,Index_jx33_km(end)) = M(end-5,Index_jx33_km(end)) - 1i/epsilon0; % for the second item of Ex

% for Ey
Index_jy11_km = getIndexOfBlkMatrix_mixed(S_km,S_bm,N,J,kappas_km,8);
Index_jy21_km = getIndexOfBlkMatrix_mixed(S_km,S_bm,N,J,kappas_km,9);
Index_jy31_km = getIndexOfBlkMatrix_mixed(S_km,S_bm,N,J,kappas_km,10);
Index_jy12_km = getIndexOfBlkMatrix_mixed(S_km,S_bm,N,J,kappas_km,11);
Index_jy22_km = getIndexOfBlkMatrix_mixed(S_km,S_bm,N,J,kappas_km,12);
Index_jy32_km = getIndexOfBlkMatrix_mixed(S_km,S_bm,N,J,kappas_km,13);
Index_jy33_km = getIndexOfBlkMatrix_mixed(S_km,S_bm,N,J,kappas_km,14);
M(end-4,Index_jy11_km(end)) = M(end-4,Index_jy11_km(end)) - 1i/epsilon0; % for the third item of Ey
M(end-4,Index_jy21_km(end)) = M(end-4,Index_jy21_km(end)) - 1i/epsilon0; % for the third item of Ey
M(end-4,Index_jy31_km(end)) = M(end-4,Index_jy31_km(end)) - 1i/epsilon0; % for the third item of Ey
M(end-4,Index_jy12_km(end)) = M(end-4,Index_jy12_km(end)) - 1i/epsilon0; % for the third item of Ey
M(end-4,Index_jy22_km(end)) = M(end-4,Index_jy22_km(end)) - 1i/epsilon0; % for the third item of Ey
M(end-4,Index_jy32_km(end)) = M(end-4,Index_jy32_km(end)) - 1i/epsilon0; % for the third item of Ey
M(end-4,Index_jy33_km(end)) = M(end-4,Index_jy33_km(end)) - 1i/epsilon0; % for the third item of Ey

% Ez
Index_jz11_km = getIndexOfBlkMatrix_mixed(S_km,S_bm,N,J,kappas_km,15);
Index_jz21_km = getIndexOfBlkMatrix_mixed(S_km,S_bm,N,J,kappas_km,16);
Index_jz31_km = getIndexOfBlkMatrix_mixed(S_km,S_bm,N,J,kappas_km,17);
Index_jz12_km = getIndexOfBlkMatrix_mixed(S_km,S_bm,N,J,kappas_km,18);
Index_jz22_km = getIndexOfBlkMatrix_mixed(S_km,S_bm,N,J,kappas_km,19);
Index_jz32_km = getIndexOfBlkMatrix_mixed(S_km,S_bm,N,J,kappas_km,20);
Index_jz13_km = getIndexOfBlkMatrix_mixed(S_km,S_bm,N,J,kappas_km,21);
Index_jz23_km = getIndexOfBlkMatrix_mixed(S_km,S_bm,N,J,kappas_km,22);
Index_jz33_km = getIndexOfBlkMatrix_mixed(S_km,S_bm,N,J,kappas_km,23);
Index_jz34_km = getIndexOfBlkMatrix_mixed(S_km,S_bm,N,J,kappas_km,24);
M(end-3,Index_jz11_km(end)) = M(end-3,Index_jz11_km(end)) - 1i/epsilon0; % for the second item of Ez
M(end-3,Index_jz21_km(end)) = M(end-3,Index_jz21_km(end)) - 1i/epsilon0; % for the second item of Ez
M(end-3,Index_jz31_km(end)) = M(end-3,Index_jz31_km(end)) - 1i/epsilon0; % for the second item of Ez
M(end-3,Index_jz12_km(end)) = M(end-3,Index_jz12_km(end)) - 1i/epsilon0; % for the second item of Ez
M(end-3,Index_jz22_km(end)) = M(end-3,Index_jz22_km(end)) - 1i/epsilon0; % for the second item of Ez
M(end-3,Index_jz32_km(end)) = M(end-3,Index_jz32_km(end)) - 1i/epsilon0; % for the second item of Ez
M(end-3,Index_jz13_km(end)) = M(end-3,Index_jz13_km(end)) - 1i/epsilon0; % for the second item of Ez
M(end-3,Index_jz23_km(end)) = M(end-3,Index_jz23_km(end)) - 1i/epsilon0; % for the second item of Ez
M(end-3,Index_jz33_km(end)) = M(end-3,Index_jz33_km(end)) - 1i/epsilon0; % for the second item of Ez
M(end-3,Index_jz34_km(end)) = M(end-3,Index_jz34_km(end)) - 1i/epsilon0; % for the second item of Ez

%% Step 7: The perturbed currents of Maxwell's equations for the bi-Maxwellian plasmas
% for Ex
Index_jx11_bm = getIndexOfBlkMatrix_mixed(S_km,S_bm,N,J,kappas_km,25);
Index_jx21_bm = getIndexOfBlkMatrix_mixed(S_km,S_bm,N,J,kappas_km,26);
Index_jx31_bm = getIndexOfBlkMatrix_mixed(S_km,S_bm,N,J,kappas_km,27);
M(end-5,end-1) = M(end-5,end-1) + c2*kz; % for the first item of Ex
M(end-5,Index_jx11_bm(end)) = M(end-5,Index_jx11_bm(end)) - 1i/epsilon0; % for the second item of Ex
M(end-5,Index_jx21_bm(end)) = M(end-5,Index_jx21_bm(end)) - 1i/epsilon0; % for the second item of Ex
M(end-5,Index_jx31_bm(end)) = M(end-5,Index_jx31_bm(end)) - 1i/epsilon0; % for the second item of Ex
M(end-5,end-8) = M(end-5,end-8) - 1i/epsilon0; % for the second item of Ex with jx

% for Ey
Index_jy11_bm = getIndexOfBlkMatrix_mixed(S_km,S_bm,N,J,kappas_km,28);
Index_jy21_bm = getIndexOfBlkMatrix_mixed(S_km,S_bm,N,J,kappas_km,29);
Index_jy31_bm = getIndexOfBlkMatrix_mixed(S_km,S_bm,N,J,kappas_km,30);
M(end-4,end-2) = M(end-4,end-2) - c2*kz; % for the first item of Ey
M(end-4,end) = M(end-4,end) + c2*kx; % for the second item of Ey
M(end-4,Index_jy11_bm(end)) = M(end-4,Index_jy11_bm(end)) - 1i/epsilon0; % for the third item of Ey
M(end-4,Index_jy21_bm(end)) = M(end-4,Index_jy21_bm(end)) - 1i/epsilon0; % for the third item of Ey
M(end-4,Index_jy31_bm(end)) = M(end-4,Index_jy31_bm(end)) - 1i/epsilon0; % for the third item of Ey
M(end-4,end-7) = M(end-4,end-7) - 1i/epsilon0; % for the third item of Ey

% for Ez
Index_jz11_bm = getIndexOfBlkMatrix_mixed(S_km,S_bm,N,J,kappas_km,31);
Index_jz21_bm = getIndexOfBlkMatrix_mixed(S_km,S_bm,N,J,kappas_km,32);
Index_jz31_bm = getIndexOfBlkMatrix_mixed(S_km,S_bm,N,J,kappas_km,33);
M(end-3,end-1) = M(end-3,end-1) - c2*kx;  % for the first item of Ez
M(end-3,Index_jz11_bm(end)) = M(end-3,Index_jz11_bm(end)) - 1i/epsilon0; % for the second item of Ez
M(end-3,Index_jz21_bm(end)) = M(end-3,Index_jz21_bm(end)) - 1i/epsilon0; % for the second item of Ez
M(end-3,Index_jz31_bm(end)) = M(end-3,Index_jz31_bm(end)) - 1i/epsilon0; % for the second item of Ez
M(end-3,end-6) = M(end-3,end-6) - 1i/epsilon0; % for the second item of Ez with jz

%% Step 6: Maxwell's equations for the perturbed quantities of Bx, By and Bz.
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
[wi,index]=sort(imag(omega),'descend');
% [wr,ind]=sort(real(omega),'descend');
w = omega(index);
w(abs(w)==min(abs(w))) = NaN+1i*NaN; % remove zero solution

% get Ex, Ey, Ez, Bx, By, and Bz
V_re = real(V(:,index));
V_im = imag(V(:,index));
dEx = V(end-5,index);
dEy = V(end-4,index);
dEz = V(end-3,index);
dBx = V(end-2,index);
dBy = V(end-1,index);
dBz = V(end,index);

end





