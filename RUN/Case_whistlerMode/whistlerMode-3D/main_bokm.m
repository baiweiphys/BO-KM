% filename: main.m
% main routine for oblique plasma wave model with a mixed distribution 
% in both of kappa-Maxweillian and bi-Maxwellian plasmas.
% When the plasma propagation angle is 90 degrees (deg=90), the program 
% initiates the module (PerpSolver4BK), to compute the dielectric tensor
% for propagation propagation in bi-kappa plasmas.
% Created by baiwei (email:baiweiphys@gmail.com, baiwei12@mail.ustc.edu.cn) 
% and modified by baiwei on Oct 17th, 2023

clear all;
clc;
close all;

%% input parameters
N = 5; 
J = 8; 
nth = 60;
deg = linspace(1e-4,90-1e-4,nth);
theta = deg*pi/180;
sp = 1; % sp=0: eig(); sp=1: sparse eigs();
B0 = 1.0e-6; % background magnetic field in z direction 
nk = 50;
kk0 = linspace(1e-4,0.3,nk);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (abs(deg-90)<1e-9)
    modules_path = '../../../PerpSolver4BK';
    addpath(modules_path);
    par = importdata('./bokm.in', ' ', 1); % read input parameters
    [S,kappas,Ts_parallel,Ts_perp,vts_parallel,vts_perp,wps,wcs,us0,rhocs,lambdaDs] = getPlasmaPrameters(par,B0);
    solver = @(kx,kz,theta) solver_bikappa(S,N,kx,kz,B0,par,sp);
else
    modules_path = '../../../kernels';
    addpath(modules_path);
    par = importdata('./bokm.in', ' ', 1); % read input parameters
    [S,index_skm,index_sbm,kappas,vts_parallel,vts_perp,Ts_parallel,Ts_perp,wps,wcs,us0,rhocs,lambdaDs] = getPlasmaPrameters(B0,par);
    if sum(index_sbm)==0
        solver = @(kx,kz,theta) solver_km(N,kx,kz,theta,B0,par,sp);
    elseif sum(index_skm)==0
        solver = @(kx,kz,theta) solver_maxwell(N,J,kx,kz,theta,B0,par,sp);
    else
        solver = @(kx,kz,theta) solver_mixed(N,J,kx,kz,theta,B0,par,sp);
    end
end
% [S,col]=size(par.data);
% if(col~=8)
%     disp('Wrong input data!');
% end
params_with_unit;

%%
tic;
for ith = 1:nth
    ith
    for ik = 1:nk
        k = kk0(ik)/rhocs(1);
        kx = k*sin(theta(ith)); % perpendicular to the magnetic B0
        kz = k*cos(theta(ith)); % parallel to the magnetic B0
        w = double(solver(kx,kz,theta(ith)));
        ww(ith,ik,:) = w(:);
    end
end

rmpath(modules_path); % Remove folders from search path
kk = kk0/rhocs(1);
kxx = k*sin(theta); % perpendicular to the magnetic B0
kzz = k*cos(theta); % parallel to the magnetic B0

runtime = toc;
disp(['Time elapsed: ', num2str(runtime/60), ' minutes']);
%%
for s=1:S
    ms(s)=par.data(s,2)*me; % mass
    ns0(s)=par.data(s,3); % desity unit: m^-3
    Ts_parallel(s) = par.data(s,4)*qe/kB; % parallel temperature, unit: eV -> K
    Ts_perp(s) = par.data(s,5)*qe/kB; % perp temperature, unit: eV -> K
end

betas_para = 2*mu0*kB.*ns0.*Ts_parallel./B0^2; % beta_para
betas_perp = 2*mu0*kB.*ns0.*Ts_perp./B0^2; % beta_perp
vA = B0/sqrt(mu0*sum(ms.*ns0)); % Alfven speed
%
disp(['beta_para = ', num2str(betas_para)]);
disp(['beta_perp = ', num2str(betas_perp)]);
disp(['Alfven speed = ', num2str(vA)]);

%% Save data
% Create a new data folder if it doesn't exist
currentPath = pwd;
foldername = 'output';
run('../../../toolBox/createDateFile(currentPath,foldername)');
fname = ['./',foldername,'/bokmData.mat'];
save(fname,'theta','kk0','kk','kk','ww','wps','wcs','lambdaDs','rhocs','N','runtime');

%% 
foldername = 'output';
fname = ['./',foldername,'/bokmData.mat'];
load(fname);
real_w = real(ww);
% real_w(real_w<0.0) = 0; 
imag_w = imag(ww);

[THETA,KK] = meshgrid(theta,kk);
KKX = KK.*sin(THETA);
KKZ= KK.*cos(THETA);

%% plot all roots
close all;
h=figure('unit','normalized','Position',[0.01 0.45 0.3 0.65],'DefaultAxesFontSize',20);
subplot(2,1,1)
NN = 100; 

for ii=36:38
    rww = reshape(real_w(:,:,ii)',3000,1);
    rkk = reshape(KK*rhocs(1),3000,1);
    rth = reshape(THETA*180/pi,3000,1);
    contourf(KK*rhocs(1),THETA*180/pi,real_w(:,:,ii)'/abs(wcs(1)),[0.1,0.7],'LineColor','none');
    hold on;
    %scatter(rkk,rth,[],rww/abs(wcs(1)),'filled','MarkerFaceAlpha',0.6);
end
shading interp
colormap('jet');
colorbar;
xlabel('k\rho_{ce}');
ylabel('\phi');
% xlim([1e-3 1e5]);
ylim([0 90]);


subplot(2,1,2);
for ii=35:38
    contourf(KK*rhocs(1),THETA*180/pi,imag_w(:,:,ii)'/abs(wcs(1)),NN,'LineColor','none');
    hold on;
end
shading interp
colormap('jet');
colorbar;
xlabel('k\rho_{ce}');
ylabel('\phi');
% xlim([1e-3 1e5]);
% ylim([0 0.7]);




%% plot the maximum value of growth rate
for rootNo = 36;
iTH = 30;
close all;
h=figure('unit','normalized','Position',[0.01 0.45 0.4 0.5],'DefaultAxesFontSize',20);
subplot(2,1,1);
plot(kk*rhocs(1),real_w(iTH,:,rootNo)/abs(wcs(1)),'k.','markersize',15);
grid on;
xlabel('k\rho_{ce}');
ylabel('\omega_r/\omega_{ce}');
xlim([0 0.3]);
ylim([0 0.7]);

subplot(2,1,2);
plot(kk*rhocs(1),imag_w(iTH,:,rootNo,:)/abs(wcs(1)),'b.-','markersize',15);
grid on;
xlabel('k\rho_{ce}');
ylabel('\gamma/\omega_{ce}');
xlim([0 0.3]);
ylim([-0.05,0.0]);

pause(1);
end
