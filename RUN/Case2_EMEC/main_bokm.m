% @Description: main routine for oblique plasma wave model with a mixed 
% distribution in both of kappa-Maxweillian and bi-Maxwellian plasmas.
% When the plasma propagation angle is 90 degrees (deg=90), the program 
% calls the module (PerpSolver4BK), to compute the dielectric tensor
% for propagation propagation in bi-kappa plasmas..
% @Filename: main_bokm.m
% @Author: Bai Wei (baiweiphys@gmail.com, baiwei12@mail.ustc.edu.cn)
% @Date: 2023-08-1
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

clear all;
clc;
close all;

%% input parameters
N = 3; 
J = 8; 
deg = 1.0e-9;
%
theta = deg*pi/180;
sp = 1; % sp=0: eig(); sp=1: sparse eigs();
B0 = 1.0e-5;  % background magnetic field in z direction 
nk = 100;
kk0 = linspace(1e-4,3.0,nk);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (abs(deg-90)<1e-8)
    modules_path = '../../PerpSolver4BK';
    addpath(modules_path);
    par = importdata('./bokm.in', ' ', 1); % read input parameters
    [S,kappas,Ts_parallel,Ts_perp,vts_parallel,vts_perp,wps,wcs,us0,rhocs,lambdaDs] = getPlasmaPrameters(par,B0);
    solver = @(kx,kz,theta) solver_bikappa(S,N,kx,kz,B0,par,sp);
else
    modules_path = '../../kernels';
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
    w = double(solver(kx,kz,theta));
    ww(ik,:) = w;
    www(ik,1,:) = w;
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
run('../../toolBox/createDateFile(currentPath,foldername)');
fname = ['./',foldername,'/bokmData.mat'];
save(fname,'kk0','kk','ww','www','wps','wcs','lambdaDs','rhocs','N','runtime');

%% Load data
foldername = 'output';
fname = ['./',foldername,'/bokmData.mat'];
load(fname);
real_w = real(ww);
imag_w = imag(ww);

%% plot all result
h=figure('unit','normalized','Position',[0.01 0.45 0.5 0.6],'DefaultAxesFontSize',20);
subplot(2,1,1);
plot(kk0,real_w(:,:)/abs(wcs(1)),'k.');
xlabel('ck/\omega_{pe}');
ylabel('\omega_r/\omega_{ce}');
grid on;
xlim([0 3]);
ylim([0 1]);

subplot(2,1,2);
plot(kk0,imag_w(:,:)/abs(wcs(1)),'k.-');
xlabel('ck/\omega_{pe}');
ylabel('\gamma/\omega_{ce}');
grid on;
xlim([0 3]);
ylim([-0.05,0.2]);

%% get the maximun value of growth rate
rootNo = 1:37;

h=figure('unit','normalized','Position',[0.01 0.45 0.5 0.6],'DefaultAxesFontSize',20);
subplot(2,1,1);
plot(kk0,real_w(:,rootNo)/abs(wcs(1)),'k.','markersize',15);
xlabel('ck/\omega_{pe}');
ylabel('\omega_r/\omega_{ce}');
grid on;
xlim([0 3]);
ylim([0 1]);


subplot(2,1,2);
plot(kk0,imag_w(:,rootNo)/abs(wcs(1)),'b.-','markersize',15);
xlabel('ck/\omega_{pe}');
ylabel('\gamma/\omega_{ce}');
grid on;
xlim([0 3]);
ylim([-0.05,0.2]);
