
clear all;

params_with_unit;
% read input parameters
par = importdata('input.in', ' ', 1); 

B0 = 1e-6;
[S,kappas,vts_parallel,vts_perp,wps,wcs,us0,rhocs,lambdaDs] = getPlasmaPrameters(par,B0); 
s = 1;
ms(s)=par.data(s,2)*me; % mass
qs(s)=par.data(s,1)*qe; % charge

vtsPerp2_c2 = 0.1^2;
Ts_perp = 0.5*ms(1)*c2*vtsPerp2_c2/qe

% 
wpe2_wce2 = 0.5^2;
ne = wpe2_wce2*epsilon0*B0^2/me