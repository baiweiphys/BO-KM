% function plot_all

load('saveData_benchmark_EMIC_KUPDAP_th30.mat');

for ii=1:3
    delta_Ex(:,ii) = wws(:,2,ii);
    delta_Ey(:,ii) = wws(:,3,ii);
    delta_Ez(:,ii) = wws(:,4,ii);
    delta_Bx(:,ii) = wws(:,5,ii);
    delta_By(:,ii) = wws(:,6,ii);
    delta_Bz(:,ii) = wws(:,7,ii);
end

kk0 = pas;

h=figure('unit','normalized','Position',[0.01 0.45 0.5 0.7],'DefaultAxesFontSize',20);
rootNo = 1:3;
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





