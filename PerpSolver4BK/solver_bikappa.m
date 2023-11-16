function w = solver_bikappa(S,N,kx,kz,B0,par,sp)
% @Description: To calculate the roots for the perpendicular propagation
% plasma wave model with a bi-kappa distribution, given the value of kx and kz.
% @Filename: solver_bikappa.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

params_with_unit;

[S,kappas,Ts_parallel,Ts_perp,vts_parallel,vts_perp,wps,wcs,us0,rhocs,lambdaDs] = getPlasmaPrameters(par,B0);
lambdas = 0.5*kx^2*rhocs.^2; % for argument of the modified bessel function

S0sn = @(s,n) S0snFunc(s,n,kappas,lambdas);
S1sn = @(s,n) S1snFunc(s,n,kappas,lambdas);
S3sn = @(s,n) S3snFunc(s,n,kappas,lambdas);
S5sn = @(s,n) S5snFunc(s,n,kappas,lambdas);

% for x-component
bx1sn_bikappa = @(s,n) biKappa_bx1sn(s,n,kappas,wps,lambdas,S1sn);
bx2sn_bikappa = @(s,n) biKappa_bx2sn(s,n,kappas,wps,lambdas,S3sn);

% for y-component
by1sn_bikappa = @(s,n) biKappa_by1sn(s,n,kappas,wps,lambdas,S3sn);
by2sn_bikappa = @(s,n) biKappa_by2sn(s,n,kappas,wps,lambdas,S5sn);

% for z-component
bz3_bikappa = biKappa_bz3(S,N,kappas,Ts_parallel,Ts_perp,wps,lambdas,S0sn);
bz3sn_bikappa = @(s,n) biKappa_bz3sn(s,n,kappas,Ts_parallel,Ts_perp,wps,lambdas,S0sn);


% Assembly matrix
%%%%%%%%%%%%%%%%%%%
% jzNo = 6 for jz
% ExyzNo = 5 for Ex
% ExyzNo = 4 for Ey
% ExyzNo = 3 for Ez
% BxyzNo = 2 for Bx
% BxyzNo = 1 for By
% BxyzNo = 0 for Bz
%%%%%%%%%%%%%%%%%%%

% Step 1
% for x-component
Mx1_bikappa = M_bikappa(S,N,wcs,bx1sn_bikappa,1,5); % No.1
Mx2_bikappa = M_bikappa(S,N,wcs,bx2sn_bikappa,2,4); % No.2

% for y-component
My1_bikappa = M_bikappa(S,N,wcs,by1sn_bikappa,3,5); % No.3
My2_bikappa = M_bikappa(S,N,wcs,by2sn_bikappa,4,4); % No.4

% for z-component
Mz_bikappa = M_bikappa(S,N,wcs,bz3sn_bikappa,5,3);  % No.5


% Step 2
MatrixColLen = size(Mx1_bikappa,2);
O = zeros(7,MatrixColLen);

M = [Mx1_bikappa;Mx2_bikappa;My1_bikappa;My2_bikappa;Mz_bikappa;O];


% for jz
M(end-6,end-3) = M(end-6,end-3) + bz3_bikappa; % for the third term of Ez

% for Ex
endIndex_jx1 = getIndexOfBlkMatrix_bikappa(S,N,1);
endIndex_jx2 = getIndexOfBlkMatrix_bikappa(S,N,2);
M(end-5,end-1) = M(end-5,end-1) + c2*kz; % for the first item of Ex
M(end-5,endIndex_jx1(end)) = M(end-5,endIndex_jx1(end)) - 1i/epsilon0; % for the second item of Ex
M(end-5,endIndex_jx2(end)) = M(end-5,endIndex_jx2(end)) - 1i/epsilon0; % for the second item of Ex

% for Ey
endIndex_jy1 = getIndexOfBlkMatrix_bikappa(S,N,3);
endIndex_jy2 = getIndexOfBlkMatrix_bikappa(S,N,4);
M(end-4,end-2) = M(end-4,end-2) - c2*kz; % for the first item of Ey
M(end-4,end) = M(end-4,end) + c2*kx; % for the second item of Ey
M(end-4,endIndex_jy1(end)) = M(end-4,endIndex_jy1(end)) - 1i/epsilon0; % for the third item of Ey
M(end-4,endIndex_jy2(end)) = M(end-4,endIndex_jy2(end)) - 1i/epsilon0; % for the third item of Ey

% Ez
endIndex_jz = getIndexOfBlkMatrix_bikappa(S,N,5);
M(end-3,end-1) = M(end-3,end-1) - c2*kx;  % for the first item of Ez
M(end-3,endIndex_jz(end)) = M(end-3,endIndex_jz(end)) - 1i/epsilon0; % for the second item of Ez
M(end-3,end-6) = M(end-3,end-6) - 1i/epsilon0; % for the second item of Ez with jz

% for Bx, By and Bz
M(end-2,end-4) = M(end-2,end-4) - kz; % for Bx
M(end-1,end-5) = M(end-1,end-5) + kz; % for the first item of By
M(end-1,end-3) = M(end-1,end-3) - kx; % for the second item of By
M(end,end-4) = M(end,end-4) + kx;     % for Bz

%% %%%%%%%% Solver
if(sp==0)
    omega = vpa(eig(M),16);
else
    omega = vpa(eigs(sparse(M),size(M,1),eps),16); 
end

[wi,ind]=sort(imag(omega),'descend');
% [wr,ind]=sort(real(omega),'descend');
w = omega(ind);
w(abs(w)==min(abs(w))) = NaN+1i*NaN; % remove zero solution


