function [Q_Short,Qp,Q_Long]=bigleaf_long_shortwave(LAT, R_inc, Ta, Tl, DOY, tm)

%Author: Giulia Vico
%Date: September 2020
%This is part of the codes used to generate the results in the following pubblication.
%Please refer to that for details on the model, its rationale and parameterization.
%Luan, X. and Vico, G. (2021), Canopy temperature and heat stress are increased by compound high air temperature
%and water stress, and reduced by irrigation ? A modeling analysis, Hydrol. Earth Syst. Sci.,
%https://doi.org/10.5194/hess-2020-549

%Key additional reference:
%Tuzet, A. et al: A coupled model of stomatal conductance, photosynthesis and transpiration, 
%Plant Cell Env, 26, 1097-1116, https://doi.org/10.1046/j.1365-3040.2003.01035.x, 2003.

%INPUTS:
%LAT latitude (in deg)
%R_inc incoming radiation (W m^(-2))
%Ta air temperature (deg C)
%Tl leaf temperature (deg C)
%DOY day of the year of interest
%tm solar time of interest (hrs)

%OUTPUTS:
%Q_short: shortwave incoming radiation (W/m2leaf)
%Qp: photosynthetically active radiation (mu mol/m2/s)
%Q_long: longwave incoming radiation (W/m2leaf) - (not used anywhere else other than to calcultae Rabs for plotting)


%FUNCTIONS CALLED: 
%CONSTANTS
%PARAMETERS
%solar_dec_azimuth




CONSTANTS;
PARAMETERS;  %it contains LAI

%--- transforms temperatures in K, if necessary
if Ta<270
    Ta=Ta+273.15;
end
if Tl<270
    Tl=Tl+273.15;
end

%---------- Zenith angle based on time of day, day of year, and site characteristics
[ZEN, AZIM, dd]=solar_dec_azimuth (DOY, LAT, tm);

% --- Compute incident radiation
R = max(0, R_inc*cos(ZEN));


%--- Compute PAR and NIR flux from total solar radiation R
PAR_top = PAR_frac*R; %[W/m2]
NIR_top = NIR_frac*R; %[W/m2]

% --- Reflection coefficients for horizontal leaves
rho_ch_PAR = (1-sqrt(1-sigm_PAR))./(1+sqrt(1-sigm_PAR));
rho_ch_NIR = (1-sqrt(1-sigm_NIR))./(1+sqrt(1-sigm_NIR));

% --- Reflection coefficients corrected for spherical leaf angle distribution
kb = 1/2/cos(ZEN); % for spherical leaf angle distribution
rho_PAR = rho_ch_PAR*2*kb./(kb+kd);
rho_NIR = rho_ch_NIR*2*kb./(kb+kd);

%--- Extinction coefficients
kestinctPAR=kb*sqrt(1-sigm_PAR);
kestinctNIR=kb*sqrt(1-sigm_NIR);

% --- Incoming shortwave components
%following Tuzet et al (2003)
Q_PAR=PAR_top*(1-rho_PAR)*(1-exp(-kestinctPAR*LAI));
Q_NIR=NIR_top*(1-rho_NIR)*(1-exp(-kestinctNIR*LAI));
Qp = q*Q_PAR;                           % photosynthetically active radiation [umol/m2/s]
Q_Short = Q_PAR + Q_NIR;                % total absorbed shortwave irradiance [W/m2_leaf]

%--- Net longwave irradiance
kestinctLW=kd;
epsil_clear = 9.2*(273.16+Ta(length(Ta)))^2/10^6;
epsil_a = epsil_clear.*(1-0.84*cloud_frac)+0.84*cloud_frac; %(Campbell and Norman, 1998, p 164)
Q_Long=(epsil_a*sigm*Ta^4-epsil*sigm*Tl^4)*(1-exp(-kestinctLW*LAI));


return

