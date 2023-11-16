function [S,kappas,Ts_parallel,Ts_perp,vts_parallel,vts_perp,wps,wcs,us0,rhocs,lambdaDs] = getPlasmaPrameters(par,B0) 
% filename: readInputPrameters_allInOne.m
% Created by Bai Wei on Jan 5th, 2022.
% Modified on Aug 16th, 2023

params_with_unit;

[S,~] = size(par.data);

for s=1:S
    qs(s)=par.data(s,1)*qe; % charge
    ms(s)=par.data(s,2)*me; % mass
    ns0(s)=par.data(s,3); % desity unit: m^-3
    Ts_parallel(s) = par.data(s,4)*qe/kB; % parallel temperature, unit: eV -> K
    Ts_perp(s) = par.data(s,5)*qe/kB; % perp temperature, unit: eV -> K
    kappas(s)=par.data(s,7);
    kappas_threshold(s)=par.data(s,8);
    As(s) = par.data(s,5)./par.data(s,4) - 1; % The s species temperature anisotropy
    
    if (kappas(s)<kappas_threshold(s))
        % thermal velocity of bi-kappa distribution
        vts_parallel(s) = sqrt(2*(1-1.5/kappas(s))*kB*Ts_parallel(s)/ms(s)); 
        vts_perp(s) = sqrt(2*kB*Ts_perp(s)/ms(s));
    else
        vts_parallel(s) = sqrt(2*kB*Ts_parallel(s)/ms(s)); % biMaxwell distribution
        vts_perp(s) = sqrt(2*kB*Ts_perp(s)/ms(s)); 
    end
    us0(s) = par.data(s,6).*sqrt(c2); % parallel drift velocity
end


lambdaDs=sqrt(epsilon0*kB*Ts_parallel./(ns0.*qs.^2)); % Debye length, Tzs
kDs=1./lambdaDs;
wps=sqrt(ns0.*qs.^2./ms/epsilon0); % plasma frequency
wcs = B0*qs./ms; % % cyclotron frequency
wce = abs(wcs(1));
rhocs = vts_perp./abs(wcs); % cyclotron radius
wps2=wps.^2;

% 
% betasz=2*mu0*kB.*ns0.*Tzs./B0^2;
% betasp=2*mu0*kB.*ns0.*Tps./B0^2;
% vA=B0/sqrt(mu0*sum(ms.*ns0)); %
% cS=sqrt(2*min(kB*Tzs)/max(ms));

