function w = solver_maxwell(N,J,kx,kz,theta,B0,par,sp)
% filename: solver_maxwell.m
% To calcuate the roots by given k for oblique plasma wave model 
% with a bi-Maxweillian distrubution.
% The script created by Bai Wei (email:baiwei12@mail.ustc.edu.cn) 
% was created on Sep 3rd, 2023，
% and revised on Sep 4th, 2023.

params_with_unit;

[S,index_skm,index_sbm,kappas,vts_parallel,vts_perp,Ts_parallel,Ts_perp,wps,wcs,us0,rhocs,lambdaDs] = getPlasmaPrameters(B0,par); 
lambdas = 0.5*kx^2*rhocs.^2; % for argument of the modified bessel function

bsnj = @(s,n,jj) maxwell_bsnj(s,n,jj,vts_parallel,Ts_parallel,Ts_perp,wcs,kz,J);
csnj = @(s,n,jj) maxwell_csnj(s,n,jj,vts_parallel,wcs,us0,kz,J);

% for x-component
bx10_maxwell = maxwell_bx10(S,N,J,bsnj,csnj,wps,lambdas);
bx11snj = @(s,n,jj) maxwell_bx11snj(s,n,jj,bsnj,csnj,wps,lambdas);
bx20_maxwell = maxwell_bx20(S,N,J,bsnj,csnj,wps,lambdas);
bx21snj = @(s,n,jj) maxwell_bx21snj(s,n,jj,bsnj,csnj,wps,lambdas);
bx30_maxwell = maxwell_bx30(S,N,J,bsnj,csnj,wps,lambdas,theta);
bx31snj = @(s,n,jj) maxwell_bx31snj(s,n,jj,bsnj,csnj,wps,wcs,lambdas,theta);


% for y-component
by10_maxwell = maxwell_by10(S,N,J,bsnj,csnj,wps,lambdas);
by11snj = @(s,n,jj) maxwell_by11snj(s,n,jj,bsnj,csnj,wps,lambdas);
by20_maxwell = maxwell_by20(S,N,J,bsnj,csnj,wps,lambdas);
by21snj = @(s,n,jj) maxwell_by21snj(s,n,jj,bsnj,csnj,wps,lambdas);
by30_maxwell = maxwell_by30(S,N,J,bsnj,csnj,wps,lambdas,theta);
by31snj = @(s,n,jj) maxwell_by31snj(s,n,jj,bsnj,csnj,wps,wcs,lambdas,theta);

% for z-component
bz10_maxwell = maxwell_bz10(S,N,J,bsnj,csnj,wps,lambdas,theta);
bz11snj = @(s,n,jj) maxwell_bz11snj(s,n,jj,bsnj,csnj,wps,wcs,lambdas,theta);
bz20_maxwell = maxwell_bz20(S,N,J,bsnj,csnj,wps,lambdas,theta);
bz21snj = @(s,n,jj) maxwell_bz21snj(s,n,jj,bsnj,csnj,wps,wcs,lambdas,theta);
bz30_maxwell = maxwell_bz30(S,N,J,bsnj,csnj,wps,lambdas,theta);
bz31snj = @(s,n,jj) maxwell_bz31snj(s,n,jj,bsnj,csnj,wps,wcs,lambdas,theta);


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
Mx11_maxwell = M_maxwell(S,N,J,csnj,bx11snj,1,5); % No.1
Mx21_maxwell = M_maxwell(S,N,J,csnj,bx21snj,2,4); % No.2
Mx31_maxwell = M_maxwell(S,N,J,csnj,bx31snj,3,3); % No.3

% for y-component
My11_maxwell = M_maxwell(S,N,J,csnj,by11snj,4,5); % No.4
My21_maxwell = M_maxwell(S,N,J,csnj,by21snj,5,4); % No.5
My31_maxwell = M_maxwell(S,N,J,csnj,by31snj,6,3); % No.6

% for z-component
Mz11_maxwell = M_maxwell(S,N,J,csnj,bz11snj,7,5); % No.7
Mz21_maxwell = M_maxwell(S,N,J,csnj,bz21snj,8,4); % No.8
Mz31_maxwell = M_maxwell(S,N,J,csnj,bz31snj,9,3); % No.9


% Step 2
MatrixColLen = size(Mx11_maxwell,2);
O = zeros(9,MatrixColLen);

M = [Mx11_maxwell;Mx21_maxwell;Mx31_maxwell;...
     My11_maxwell;My21_maxwell;My31_maxwell;...
     Mz11_maxwell;Mz21_maxwell;Mz31_maxwell;O];


% for jx
M(end-8,end-5) = M(end-8,end-5) + bx10_maxwell; % for the first term of Ex
M(end-8,end-4) = M(end-8,end-4) + bx20_maxwell; % for the second term of Ey
M(end-8,end-3) = M(end-8,end-3) + bx30_maxwell; % for the third term of Ez

% for jy
M(end-7,end-5) = M(end-7,end-5) + by10_maxwell; % for the first term of Ex
M(end-7,end-4) = M(end-7,end-4) + by20_maxwell; % for the second term of Ey
M(end-7,end-3) = M(end-7,end-3) + by30_maxwell; % for the third term of Ez

% for jz
M(end-6,end-5) = M(end-6,end-5) + bz10_maxwell; % for the first term of Ex
M(end-6,end-4) = M(end-6,end-4) + bz20_maxwell; % for the second term of Ey
M(end-6,end-3) = M(end-6,end-3) + bz30_maxwell; % for the third term of Ez

% for Ex
endIndex_jx11 = getIndexOfBlkMatrix_maxwell(S,N,J,1);
endIndex_jx21 = getIndexOfBlkMatrix_maxwell(S,N,J,2);
endIndex_jx31 = getIndexOfBlkMatrix_maxwell(S,N,J,3);
M(end-5,end-1) = M(end-5,end-1) + c2*kz; % for the first item of Ex
M(end-5,endIndex_jx11(end)) = M(end-5,endIndex_jx11(end)) - 1i/epsilon0; % for the second item of Ex
M(end-5,endIndex_jx21(end)) = M(end-5,endIndex_jx21(end)) - 1i/epsilon0; % for the second item of Ex
M(end-5,endIndex_jx31(end)) = M(end-5,endIndex_jx31(end)) - 1i/epsilon0; % for the second item of Ex
M(end-5,end-8) = M(end-5,end-8) - 1i/epsilon0; % for the second item of Ex with jx


% for Ey
endIndex_jy11 = getIndexOfBlkMatrix_maxwell(S,N,J,4);
endIndex_jy21 = getIndexOfBlkMatrix_maxwell(S,N,J,5);
endIndex_jy31 = getIndexOfBlkMatrix_maxwell(S,N,J,6);
M(end-4,end-2) = M(end-4,end-2) - c2*kz; % for the first item of Ey
M(end-4,end) = M(end-4,end) + c2*kx; % for the second item of Ey
M(end-4,endIndex_jy11(end)) = M(end-4,endIndex_jy11(end)) - 1i/epsilon0; % for the third item of Ey
M(end-4,endIndex_jy21(end)) = M(end-4,endIndex_jy21(end)) - 1i/epsilon0; % for the third item of Ey
M(end-4,endIndex_jy31(end)) = M(end-4,endIndex_jy31(end)) - 1i/epsilon0; % for the third item of Ey
M(end-4,end-7) = M(end-4,end-7) - 1i/epsilon0; % for the third item of Ey

% Ez
endIndex_jz11 = getIndexOfBlkMatrix_maxwell(S,N,J,7);
endIndex_jz21 = getIndexOfBlkMatrix_maxwell(S,N,J,8);
endIndex_jz31 = getIndexOfBlkMatrix_maxwell(S,N,J,9);
M(end-3,end-1) = M(end-3,end-1) - c2*kx;  % for the first item of Ez
M(end-3,endIndex_jz11(end)) = M(end-3,endIndex_jz11(end)) - 1i/epsilon0; % for the second item of Ez
M(end-3,endIndex_jz21(end)) = M(end-3,endIndex_jz21(end)) - 1i/epsilon0; % for the second item of Ez
M(end-3,endIndex_jz31(end)) = M(end-3,endIndex_jz31(end)) - 1i/epsilon0; % for the second item of Ez
M(end-3,end-6) = M(end-3,end-6) - 1i/epsilon0; % for the second item of Ez with jz

% for Bx, By and Bz
M(end-2,end-4) = M(end-2,end-4) - kz; % for Bx
M(end-1,end-5) = M(end-1,end-5) + kz; % for the first item of By
M(end-1,end-3) = M(end-1,end-3) - kx; % for the second item of By
M(end,end-4) = M(end,end-4) + kx;     % for Bz

%% %%%%%%%% Solver
if(sp==0)
    [V,D] = eig(M);
else
    [V,D] = eigs(sparse(M),size(M,1),eps); 
end
omega = diag(D);
[wi,ind]=sort(imag(omega),'descend');
% [wr,ind]=sort(real(omega),'descend');
w = omega(ind);
w(abs(w)==min(abs(w))) = NaN+1i*NaN; % remove zero solution


