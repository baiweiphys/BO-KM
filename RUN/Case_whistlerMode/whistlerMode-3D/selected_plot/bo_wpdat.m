% 18-10-19 17:56 Hua-sheng XIE, huashengxie@gmail.com, CCF-ENN, China
% Initial data for run bo_plot_select

% Search the most close dispersion surfaces to these data.
% Initial data for find the corresponding dispersion surfaces.
% Please use bo_plot_all.m to visualize all the solutions, and then
% modify here the initial point of which mode(s) you want plot/store.

% wpdat(:,1) is pa; wpdat(:,2) is pb for 2D scan and arbitrary for 1D scan;
% wpdat(:,3) is Re or Im(omega)


wpdat=[%25,0.2,0.435i;
       %25,0.25,-0.0089i;
       25,0.25,-0.02249i;
  ];

jselc=0; % jselc=1, alway plot the most unstable dispersion surface
% wpdat=[];

%%
% new add
wwn = ww/abs(wcs(1));
jselc=0; % jselc=1, alway plot the most unstable dispersion surface
pa = theta*180/pi;
pb = kk*rhocs(1);
npa = 60;
npb = 50;
nw = 880;
ipa = 1;
ipb = 2;
iem = 3;

%%
run ./bo_plot_select;
% subplot(122);xlim([0,2]);ylim([-1e-3,2e-3]);subplot(121);xlim([0,2]);
% subplot(122);ylim([-0.5,0.5]);
