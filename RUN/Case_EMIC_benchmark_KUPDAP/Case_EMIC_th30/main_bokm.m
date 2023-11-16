% filename: main.m
% main routine for oblique plasma wave model with a mixed distribution 
% in both of kappa-Maxweillian and bi-Maxwellian plasmas.
% When the plasma propagation angle is 90 degrees (deg=90), the program 
% call the module (PerpSolver4BK), to compute the dielectric tensor
% for propagation propagation in bi-kappa plasmas.
% Created by baiwei (email:baiweiphys@gmail.com, baiwei12@mail.ustc.edu.cn) 
% and modified by baiwei on Oct 17th, 2023

clear all;
clc;
close all;

%% input parameters
N = 4; 
J = 8; 
deg = 30.0;
%
theta = deg*pi/180;
sp = 1; % sp=0: eig(); sp=1: sparse eigs();
% B0 = 140.0e-9; % background magnetic field in z direction 
B0 = 140.0e-9;
nk = 10;
kk0 = linspace(1e-5,0.08,nk);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (abs(deg-90)<1e-3)
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
for ik = 1:nk
    ik
    k = kk0(ik)/rhocs(1);
    kx = k*sin(theta); % perpendicular to the magnetic B0
    kz = k*cos(theta); % parallel to the magnetic B0
    [w,dEx,dEy,dEz,dBx,dBy,dBz] = solver(kx,kz,theta);
    ww(ik,:) = w;
    www(ik,1,:) = w;
    www(ik,2,:) = dEx;
    www(ik,3,:) = dEy;
    www(ik,4,:) = dEz;
    www(ik,5,:) = dBx;
    www(ik,6,:) = dBy;
    www(ik,7,:) = dBz;
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
save(fname,'kk0','kk','ww','www','wps','wcs','lambdaDs','rhocs','N','runtime');

%% 
foldername = 'output';
fname = ['./',foldername,'/bokmData.mat'];
load(fname);
real_w = real(ww);
imag_w = imag(ww);
for ik=1:length(kk0)
    delta_Ex(ik,:) = www(ik,2,:);
    delta_Ey(ik,:) = www(ik,3,:);
    delta_Ez(ik,:) = www(ik,4,:);
    delta_Bx(ik,:) = www(ik,5,:);
    delta_By(ik,:) = www(ik,6,:);
    delta_Bz(ik,:) = www(ik,7,:);
end

%% plot all result
h=figure('unit','normalized','Position',[0.01 0.45 0.5 0.6],'DefaultAxesFontSize',20);
subplot(2,1,1);
plot(kk0,real_w(:,:)/abs(wcs(1)),'k.','markersize',15);
ylabel('\omega_r/\omega_{ce}');
xlabel('k\rho_{ce}');
% xlim([1e-3 1e5]);
ylim([0 1]);
grid on;
hold on;

subplot(2,1,2);
plot(kk0,imag_w(:,:)/abs(wcs(1)),'b.-','markersize',15);
ylabel('\gamma/\omega_{ce}');
xlabel('k\rho_{ce}');
grid on;
% xlim([0,0.37]);
ylim([-0.1,0.1]);


%% plot the maximum value of growth rate
h=figure('unit','normalized','Position',[0.01 0.45 0.5 0.6],'DefaultAxesFontSize',20);
rootNo = 1:10;
subplot(2,1,1);
plot(kk0,real_w(:,rootNo)/abs(wcs(1)),'k.','markersize',5);
ylabel('\omega_r/\omega_{ce}');
xlabel('k\rho_{ce}');
% xlim([1e-3 1e5]);
ylim([0 1]);
grid on;

% subplot(2,1,2);
figure;
% h=figure('unit','normalized','Position',[0.01 0.45 0.7 0.45],'DefaultAxesFontSize',20);
plot(kk0,imag_w(:,rootNo)/abs(wcs(1)),'b.-','markersize',5);
ylabel('\gamma/\omega_{ce}');
xlabel('k\rho_{ce}');
grid on;
% xlim([0,0.37]);
ylim([-0.05 0.03]);

% pause(0.2);
% title(['rootNo.=',num2str(rootNo)]);

%% plot dEx, dEy, dEz, dBx, dBy, and dBz
h=figure('unit','normalized','Position',[0.01 0.45 0.5 0.7],'DefaultAxesFontSize',20);
subplot(3,4,1);
plot(kk0,real(delta_Ex(:,rootNo)),'b.','markersize',5);
ylabel('real(\delta E_x)');
xlabel('k\rho_{ce}');
grid on;
%
subplot(3,4,2);
plot(kk0,imag(delta_Ex(:,rootNo)),'b.','markersize',5);
ylabel('imag(\delta E_x)');
xlabel('k\rho_{ce}');
grid on;

%
subplot(3,4,3);
plot(kk0,real(delta_Ey(:,rootNo)),'b.','markersize',5);
ylabel('real(\delta E_y)');
xlabel('k\rho_{ce}');
grid on;

subplot(3,4,4);
plot(kk0,imag(delta_Ey(:,rootNo)),'b.','markersize',5);
ylabel('imag(\delta E_y)');
xlabel('k\rho_{ce}');
grid on;

subplot(3,4,5);
plot(kk0,real(delta_Ez(:,rootNo)),'b.','markersize',5);
ylabel('real(\delta E_z)');
xlabel('k\rho_{ce}');
grid on;
%
subplot(3,4,6);
plot(kk0,imag(delta_Ez(:,rootNo)),'b.','markersize',5);
ylabel('imag(\delta E_z)');
xlabel('k\rho_{ce}');
grid on;


%
subplot(3,4,7);
plot(kk0,real(delta_Bx(:,rootNo)),'b.','markersize',5);
ylabel('real(\delta B_x)');
xlabel('k\rho_{ce}');
grid on;
%
subplot(3,4,8);
plot(kk0,imag(delta_Bx(:,rootNo)),'b.','markersize',5);
ylabel('imag(\delta B_x)');
xlabel('k\rho_{ce}');
grid on;


%
subplot(3,4,9);
plot(kk0,real(delta_By(:,rootNo)),'b.','markersize',5);
ylabel('real(\delta B_y)');
xlabel('k\rho_{ce}');
grid on;
%
subplot(3,4,10);
plot(kk0,imag(delta_By(:,rootNo)),'b.','markersize',5);
ylabel('real(\delta B_y)');
xlabel('k\rho_{ce}');
grid on;

%
subplot(3,4,11);
plot(kk0,real(delta_Bz(:,rootNo)),'b.','markersize',5);
ylabel('real(\delta B_z)');
xlabel('k\rho_{ce}');
grid on;
%
subplot(3,4,12);
plot(kk0,imag(delta_Bz(:,rootNo)),'b.','markersize',5);
ylabel('imag(\delta B_z)');
xlabel('k\rho_{ce}');
grid on;





