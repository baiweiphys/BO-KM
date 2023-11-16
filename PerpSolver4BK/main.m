% @Description: Driver routine to calculate the solution for the bi-kappa 
% plasma dispersion relation in a perpendicular propagation plasma wave 
% model with a bi-kappa distribution.
% @Filename: main.m
% @Date: 2021-10-01
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

clear all;
clc;
close all;

params_with_unit;

% read input parameters
par = importdata('input.in', ' ', 1); 
[S,col]=size(par.data);
if(col~=8)
    disp('Wrong input data!');
end

B0 = 5.5e-2;  % background magnetic field in z direction 
[S,kappas,Ts_parallel,Ts_perp,vts_parallel,vts_perp,wps,wcs,us0,rhocs,lambdaDs] = getPlasmaPrameters(par,B0);
%%
%%%%% input parameters
N = 5; 
deg = 90.0;
theta = deg*pi/180;
sp = 1; % sp=0: eig(); sp=1: sparse eigs();

kk=[]; 
kxx=[]; 
kzz=[]; 
ww=[];

nk = 100;
kk0 = linspace(1e-3,10,nk);

tic;
for ik = 1:nk  %kj=0.01:0.01:1.0% 
    ik
    % k = kj/rhocs(1);
    k = kk0(ik)/rhocs(1);
    kx = k*sin(theta); % the component of the k perpendicular to the magnetic B0
    kz = k*cos(theta); % the component of the k parallel to the magnetic B0
    w = double(solver_bikappa(S,N,kx,kz,B0,par,sp));
    ww(ik,:) = w;
    www(ik,1,:) = w;
end
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
save('saveData_bikappa.mat','kk0','kk','ww','www','wps','wcs','lambdaDs','rhocs','N')

%% 
load('saveData_bikappa.mat');
real_w = real(ww);
imag_w = imag(ww);

%% plot
h=figure('unit','normalized','Position',[0.01 0.45 0.5 0.6],'DefaultAxesFontSize',20);
% hold on;
plot(kk*rhocs(1),real_w(:,2:end)/abs(wcs(1)),'r.','markersize',8);
ylabel('\omega_r/\omega_{ce}');
xlabel('k\rho_{ce}');
grid on;
xlim([0 10]);
ylim([0 10]);



h=figure('unit','normalized','Position',[0.01 0.45 0.5 0.6],'DefaultAxesFontSize',20);
% hold on
plot(kk*rhocs(1),imag_w(:,2:end)/abs(wcs(1)),'r.','markersize',15);
ylabel('\gamma/\omega_{ce}');
xlabel('k\rho_{ce}');
grid on;
xlim([0 10]);
ylim([0.0,1e-13]);

%%
hold on;
plot(kk*rhocs(1),real_w(:,2:end)/abs(wcs(1)),'r^','markersize',5);
% ylabel('\omega_r/\omega_{ce}');
% xlabel('k\rho_{ce}');
% grid on;
% % xlim([0 5]);
% % ylim([0 6]);
