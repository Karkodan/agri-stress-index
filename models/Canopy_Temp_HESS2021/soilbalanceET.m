function [snew,ds,Irrinput,LQlosses]=soilbalanceET(s,ETcanopy,rain,integr_step,coeff,sfc,s1,Zr,n,Ksat,bsoil,stilde,shat)

%Authors: Xiangyu Luan, Giulia Vico 
%Date: October 2020
%This is part of the codes used to generate the results in the following pubblication. 
%Please refer to that for details on the model, its rationale and parameterization.
%Luan, X. and Vico, G. (2021), Canopy temperature and heat stress are increased by compound high air temperature 
%and water stress, and reduced by irrigation ? A modeling analysis, Hydrol. Earth Syst. Sci., 
%https://doi.org/10.5194/hess-2020-549 


%INPUTS: 
% s soil moisture at previous step
% ET canopy canopy transpiration rate, mol H20/m2ground/s, transformed below into m/s
% rain rainfall input  mm (0 if no rain occurs)
% integr_step  hrs
% coeff upscales the instantaneous ET (maximum for the day) to daily totals
% sfc soil field capacity (UNUSED, beacuse of the simplified description of runoff and deep percolation)
% s1 soil moisture threshold above which losses are instantaneous (used if integr_step>8hrs)
% Zr active rooting depth (m)
% n soil porosity
% Ksat hydraulic conductivity at saturation (m/s) - used for deep percolation, if integr_step<8hrs
% b soil exponent of the soil moisture hydraulic conductivity function - used for deep percolation, if integr_step<8hrs
% stilde: soil moisture at which an irrigation application is triggered
% shat: soil moisture target level for the irrigation application (n Zr (shat-stilde) is the depth of applied water via irrigation

% OUTPUTS:
% snew: updated soil moisture value
% ds: change in soil moisture
% Irrinput: input via irrigation (m over the period)
% LQlosses: losses via runoff and deep percolation (m over the period)


%FUNCTION CALLED:
%none


%evapotranspiration losses
ET=ETcanopy*18/10^6;  %to go from mol H20/m2/s to m^3/m2/s, i.e. m/s directly
%use directly ET (and density of liquid water) to determine the change in soil moisture!

% -------------irrigation----------
if s<stilde
    Irr=shat-stilde; %expressed as a change in s, not a water depth
else
    Irr=0;
end
Irrinput=Irr*n*Zr; %m over the period considered

%--------------------------------------------------------
if integr_step<8 %hrs, i.e., if working at the subdaily time scale 
    %LQ=Ksat/(exp(beta*(1-sfc)-1))*(exp(beta*(s-sfc)-1)); %LQ is leakage, this equation is based on Laio et al.(2001), DOI: 10.1016/s0309-1708(01)00005-7
    LQ=Ksat*s^(2*bsoil+3);  %more common appraoch to describe losses via runoff and deep percolation - rate, in m/s
    ds=(rain/(n*Zr*1000))+Irr-(ET*coeff/(n*Zr)+LQ/(n*Zr))*integr_step*3600; 
    snew=max(min(s+ds,1),0);
    LQlosses=LQ*integr_step*3600+max(s+ds-1,0)*n*Zr; %m over the period considered
else %at the diurnal scale, 
    %instantaneous loss via LQ whenever soil moisture above a threshold s1
    ds=(rain/(n*Zr*1000))+Irr-(ET*coeff/(n*Zr))*integr_step*3600; 
    snew=max(min(s+ds,s1),0);
    LQlosses=max(s+ds-s1,0)*n*Zr; %m over the period considered
end
   
return
