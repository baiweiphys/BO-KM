% 18-10-19 17:56 Hua-sheng XIE, huashengxie@gmail.com, CCF-ENN, China
% Initial data for run bo_plot_select

% Search the most close dispersion surfaces to these data.
% Initial data for find the corresponding dispersion surfaces.
% Please use bo_plot_all.m to visualize all the solutions, and then
% modify here the initial point of which mode(s) you want plot/store.

% wpdat(:,1) is pa; wpdat(:,2) is pb for 2D scan and arbitrary for 1D scan;
% wpdat(:,3) is Re or Im(omega)


wpdat=[0.0137457,0,0.00252279i;
       0.0121297,0,-0.00066583i;
       0.0161696,0,-0.00110061i;
  ];


jselc=0; % jselc=1, alway plot the most unstable dispersion surface
% wpdat=[];

%%
% new add

npa = 100;
npb = 7;
nw = 3327;
ipa = 1;
ipb = 1;
iloga = 0;
jsetplot = 1;
iem = 3;

wwn = www/abs(wcs(1));
pb = 0;
pa = kk0;


%%
run ./bo_plot_select;
% subplot(122);xlim([0,2]);ylim([-1e-3,2e-3]);subplot(121);xlim([0,2]);
% subplot(122);ylim([-0.5,0.5]);
