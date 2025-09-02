function [Tlnew]=Tleaf_layers(pa,gh,Tair,D,gv,es,Qshort,csi)

%Author: Xiangyu Luan and Giulia Vico 
%Date: September 2020
%This is part of the codes used to generate the results in the following pubblication. 
%Please refer to that for details on the model, its rationale and parameterization.
%Luan, X. and Vico, G. (2021), Canopy temperature and heat stress are increased by compound high air temperature 
%and water stress, and reduced by irrigation ? A modeling analysis, Hydrol. Earth Syst. Sci., 
%https://doi.org/10.5194/hess-2020-549 

%INPUTS: 
% pa is atmospheric pressure, unit kPa
% gh total canopy conductance to heat (per unit ground area), mol/m2/s
% Tair air temperature (deg C)
% D water vapour deficit, unit kPa
% gv total canopy conductance to water vapor (per unit ground area), mol/m2/s
% Rabs total ABSORBED radiation (W/m2)
% csi amount of cumulated leaves above the layer of interest (for a big leaf model, leaf area index) 

%OUTPUTS:
%Tlnew leaf temperature (deg C) as determined by the leaf energy balance

%FUNCTIONS CALLED: 
%CONSTANTS
%PARAMETERS


CONSTANTS;
PARAMETERS;


delta=17.502*240.97*es/(240.97+Tair).^2;
ss=delta/pa; % ss slope of the vapor pressure vs. temperature curve	, unit Pa/K

epsaclear=9.2*10^-6.*(Tair+273.15).^2; % clear sky emissivity epsilon, based on Swinbank's formula (Campbell and Norman, p 164).
epsa=epsaclear*(1-0.84*cloud_frac)+0.84*cloud_frac; %Campbell and Norman 


if length(csi)>2 %i.e., if truly layered canopy is considered
    numer= Qshort+(epsa*sigm*(Tair+273.15).^4-epsil*sigm*(Tair+273.15).^4)*kd*exp(-kd*csi)-lambd*gv*D/P; %D is in kPa!
    denomin= Cp*gh+lambd*gv*ss+4*sigm*epsil*(Tair+273.15).^3*kd*exp(-kd*csi);
else %for big leaf model
    numer= Qshort+(epsa*sigm*(Tair+273.15).^4-epsil*sigm*(Tair+273.15).^4)*(1-exp(-kd*LAI))-lambd*gv*D/P; %D is in kPa!
    denomin= Cp*gh+lambd*gv*ss+4*sigm*epsil*(Tair+273.15).^3*(1-exp(-kd*LAI));
end

Tlnew=Tair+numer/denomin;


return
