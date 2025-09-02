%Authors: Gabrel G Katul, Giulia Vico 
%Date: December 2020
%Part of the codes used to generate the results in the following pubblication. 
%Please refer to that for details on the model, its rationale and parameterization.
%Luan, X. and Vico, G. (2021), Canopy temperature and heat stress are increased by compound high air temperature 
%and water stress, and reduced by irrigation ? A modeling analysis, Hydrol. Earth Syst. Sci., 
%https://doi.org/10.5194/hess-2020-549 


%------- Conversion factors and constants

%--- Mass transfer
aa = 1.65;              % conversion from CO2 stomatal conductance to H2O stomatal conductance (Bonan 2008)
bb = 1.37;              % conversion from CO2 boundary layer conductance to H2O boundary layer conductance (Bonan 2008)
k_C = 1.4*0.110;        % coefficient for boundary layer conductance to CO2 (Campbell and Norman, Tab. 7.6)
k_W = 1.4*0.147;        % coefficient for boundary layer conductance to H2O (Campbell and Norman, Tab. 7.6) = k_C*bb
k_H = 1.4*0.135;        % coefficient for boundary layer conductance to heat (Campbell and Norman, Tab. 7.6)
gg = 9.81;              % gravity [m/s2]
kvonkarman=0.41;        %von Karman constant

%--- Air and water vapor
PM_air = 28.8;          % molecular weight of air [g/mol]
PM_w = 18;              % molecular weight of water [g/mol]
rho_w = 10^3;           % liquid water density [kg/m3]
CF = 1.15*1000/PM_air;  % air molar density in standard conditions [mol/m3]
rho_molgas=CF;%44.6;        % molecular density of a gas, under standard conditions (mol/m3)
P = 101.3;              % atmospheric pressure [kPa]
Cp = 29.3;              % heat capacity of air (or molar speacifi heat of air) [J/mol/degC]
lambdKg=2.454*10^6;       %latent heat of vaporization of water (at 20 deg C); J/kg
lambd=lambdKg*18/10^3;  %latent heat of vaporization of water (at 20 deg C); J/mol, to match units
% constants for the partial pressure of water at saturation esa = a * exp((b.*T)./(T + c)) [kPa, deg C]
a = 0.611;              % [kPa]
b = 17.502;             % [unitless]
c = 240.97;             % [deg C]
Rgas=8.314;             %[J/mol/K] universal constant of gasses

%--- Energy balance
sigm = 5.67*10^-8;     % Stefan-Boltzman constant [W/m2/K4] %do not use teh symbol sigma, because it is a MatLab function 
q = 4.6;                % conversion for PAR from J to umol (Campbell and Norman, p. 149) [umol quanta/J]
PAR_frac = 0.45;        % fraction of total radiation in the PAR wavebands (Campbell and Norman, p. 161)
NIR_frac = 1-PAR_frac;  % fraction of total radiation in the NIR wavebands
epsil = 0.97;           % long-wave LEAF emissivity (Campbell and Norman, Table 11.3)
epssoil=0.94;           % long-wave soil emissivity (Campbell and Norman, Table 11.3)
rho_d_PAR = 0.057;      % reflection coefficient in PAR bands for diffuse radiation (Leuning 1995, p. 1197)
rho_d_NIR = 0.389;      % reflection coefficient in NIR bands for diffuse radiation (Leuning 1995, p. 1197)
sigm_PAR = 0.2;         % scattering coefficient in PAR bands (Leuning 1995, p. 1197)
sigm_NIR = 0.8;         % scattering coefficient in NIR bands (Leuning 1995, p. 1197)
tau_a = 1;              % air transmissivity (tau_a>0.7 for clear sky)
a_S = 0.5;              % mean short-wave leaf absorbivity (Campbell and Norman, Table 11.4)
a_PAR = 0.85;           % PAR leaf absorbivity (Bonan)
a_NIR = 0.15;           % NIR leaf absorbivity (Bonan)
a_L = epsil;            % long-wave leaf absorbivity
kd = 0.8;               % extintion coefficient for diffuse and thermal radiation assuming spherical leaf angle distribution (Goundriaan and van Laar 1994, p 103)


