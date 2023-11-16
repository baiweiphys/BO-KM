% 18-10-19 17:56 Hua-sheng XIE, huashengxie@gmail.com, CCF-ENN, China
% Initial data for run bo_plot_select

% Search the most close dispersion surfaces to these data.
% Initial data for find the corresponding dispersion surfaces.
% Please use bo_plot_all.m to visualize all the solutions, and then
% modify here the initial point of which mode(s) you want plot/store.

% wpdat(:,1) is pa; wpdat(:,2) is pb for 2D scan and arbitrary for 1D scan;
% wpdat(:,3) is Re or Im(omega)

wpdat=[1.0e-4, 0.00647384, 0.0382281; % L mode for O+ band
       1.0e-4, 0.0194015,  0.19781;   % L mode for He+ band
       1.0e-4, 0.0194015,  0.411745;  % L mode for proton band
       1.0e-4, 0.0484888,  0.734527;  % R mode
  ];

% jselc=1, alway plot the most unstable dispersion surface
jselc = 0; 
% wpdat=[];

%%
% new add
wwn = ww/abs(wcs(1));
pa = theta*180/pi;
pb = kk*rhocs(1);
npa = 60;
npb = 100;
nw = 1137;
ipa = 1;
ipb = 2;
iem = 3;

%%
run ./bo_plot_select;
% subplot(122);xlim([0,2]);ylim([-1e-3,2e-3]);subplot(121);xlim([0,2]);
% subplot(122);ylim([-0.5,0.5]);
