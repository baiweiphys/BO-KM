% plot S0, S1, S3, and S5 functions


clear all;
clc;
close all;

num = 100;
lambda = linspace(1e-4,1e-1,num);

s = 1;
n = 1;
kappas = 5;


S0sn = @(lambdas) S0snFunc(s,n,kappas,lambdas);
S1sn = @(lambdas) S1snFunc(s,n,kappas,lambdas);
S3sn = @(lambdas) S3snFunc(s,n,kappas,lambdas);
S5sn = @(lambdas) S5snFunc(s,n,kappas,lambdas);

S0 = [];
S1 = [];
S3 = [];
S5 = [];
for i=1:num
    S0 = [S0, S0sn(lambda(i))];
    S1 = [S1, S1sn(lambda(i))];
    S3 = [S3, S3sn(lambda(i))];
    S5 = [S5, S5sn(lambda(i))];
end


h=figure('unit','normalized','Position',[0.01 0.45 0.5 0.6],'DefaultAxesFontSize',20);
plot(lambda,real(S0),'LineWidth',3);
hold on;
plot(lambda,real(S1),'LineWidth',3);
plot(lambda,real(S3),'LineWidth',3);
plot(lambda,real(S5),'LineWidth',3);
legend('S_{0}','S_{1}','S_{3}','S_{5}');
% legend('boxoff');
xlabel('\lambda');
ylabel('Re(S_j)');



