function [gsrp,gsr,gp,gpalt]=soil2leaf_conductance(psil,psis,psissat,Ksat,bsoil,Zr,RAI)

%Author: Giulia Vico
%Date: September 2020
%Part of the codes used to generate the results in the following pubblication.
%Please refer to that for details on the model, its rationale and parameterization.
%Luan, X. and Vico, G. (2021), Canopy temperature and heat stress are increased by compound high air temperature
%and water stress, and reduced by irrigation ? A modeling analysis, Hydrol. Earth Syst. Sci.,
%https://doi.org/10.5194/hess-2020-549

%INPUT:
%psil: leaf water potential (MPa, with sign)
%psis: soil water potential (MPa, with sign)
%psissat: soil water potential at saturation (MPa, with sign)
%Ksat: soil hydraulic conductivity at saturation (m/s)
%bsoil: exponent of the soil water retention curve (adimensional)
%Zr: active rooting depth (m)
%LAI: leaf area index (m2leaf/m2ground)

%OUTPUT:
%gsrp: soil to root to plant conductance - m/(s MPa) i.e. m^3H20/(m^2 ground s MPa)
%gsr:  soil to root conductance
%gp: root to leaf conductance, based on Manzoni et al 2011 Advances in Water Resources
%gpalt: root to leaf conductance, alternative calculation, based on Vico and Porporato 2008 Plant and Soil

%FUNCTIONS CALLED:
%PARAMETERS
%CONSTANTS


%NOTE:
%Estimates of conductances along the SPAC are rather uncertain.
%The point of departure to determine gpmax is the leaf specific conductivity, which represents the
%conductivity normalized by the leaf area; this is different from K_sap, which is normalized by the
%sapwood area (and, as such, requires the sapwood area index - SAI - to be transformed into a per unit ground
%area quantity).
%To estimate gpmax, multiply the leaf specific conductance (in kg/(s m MPa) for LAI (m^2leaf/m^2ground),
%divide by plant height (m) and divide by water density (1000 kg/m^3). Typical values (from Kocachinar and Sage 2004)
%are: (0.5 to 6)*10^(-4) kg/(s m MPa); they become (2-24)*10^(-7) %m/(s MPa), considering a LAI of 2 and a plant height of 0.5 m
%See Traits_herbaceous_conductivities.xls for more values
%In Vico and Porporato, 2008, gpmax=1.3*10^(-7) m /(s MPa), i.e., at the lowest end of the spectrum
%NB: conductance: multiply per difference in water potential;
%conductivity: multiply by gradient of water potential (i.e., difference in water potential per unit length)

PARAMETERS;
CONSTANTS;


%SOIL TO ROOT CONDUCTANCE
Ks=Ksat*(psissat/psis)^(2*(1+1/bsoil)); %soil hydraulic conductivity at psis (m/s)
%Vico and Porporato 2008
%gsr=10^6*Ks*sqrt(RAI)/(pi*gg*rho_w*Zr); %m/(s MPa) ----- the 10^6 is to go from m/(s Pa) to m/(s MPa) in gsr
%ALTERNATIVE approach in Manzoni 2011 AWR
lsr=sqrt(dr*Zr/(RAI)); %m --- this formulation has the disadvantage that gsrp declines slower than gs, so that psil increases under certain conditions, even if psis goes down...
%if using Daly et al 2004 Journal of hydrometeorology, Part 1)
%lsr=sqrt(dr*Zr/(RAI*(psis/psissat)^(-(-Edosexponent)/bsoil))); %m
gsr=10^6*Ks/(rho_w*gg*lsr); %m/(s MPa) ----- the 10^6 is to go from m/(s Pa) to m/(s MPa) in gsr


%ROOT TO LEAF CONDUCTANCE
%Vico and Porporato 2008
gp=gpmax*exp(-(-psil/dc)^bc);
gpalt=gp;
%ALTERNATIVE approach based on Manzoni et al 2011 Advances in Water Resources
gp=gpmax/(1+(psil/psi50)^acav);



% SOIL TO LEAF CONDUCTANCE, as series
gsrp=gsr*gp/(gsr+gp); %m/(s MPa)

return