


modules_path = '../../../kernels';
addpath(modules_path);

params_with_unit;

% par = importdata('./bokm.in', ' ', 1); % read input parameters
% qs(s)=par.data(s,1)*qe; % charge
% ms(s)=par.data(s,2)*me; % mass
n0 = 6.503e6;
B0 = 173e-9;
betas_para = [8.69e-6,9.22e-5,0.000112,0.00251,0.0132,0.108,8.69e-5,8.69e-5];

Ts_parallel = betas_para*B0^2/(2*mu0*kB.*n0);


