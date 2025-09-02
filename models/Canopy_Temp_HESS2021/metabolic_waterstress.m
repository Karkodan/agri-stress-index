function [Vcmax_ws,J_ws,Rd_ws,gamma_star,k1_ws,k2_ws]=metabolic_waterstress(Vcmax25,Jmax25,psil,Qp,Tl)

%Authors: Giulia Vico 
%Date: September 2020
%This is part of the codes used to generate the results in the following pubblication. 
%Please refer to that for details on the model, its rationale and parameterization.
%Luan, X. and Vico, G. (2021), Canopy temperature and heat stress are increased by compound high air temperature 
%and water stress, and reduced by irrigation ? A modeling analysis, Hydrol. Earth Syst. Sci., 
%https://doi.org/10.5194/hess-2020-549 

%INPUT: 
% Vcmax25  maximum carboxylation rate at 25 deg C and well watered cond (mumol/m2leaf/s)
% Jmax25   light saturated rate of electron transport at 25 deg C and well watered cond (mumol/m2leaf/s)
% psil     leaf water potential, with sign (MPa)
% Qp       photosynthetically active radiation (mumol m-2 s-1)
% Tl       leaf or canopy temperature (deg C)

%OUTPUT:
%the subscript ws refers to water stressed conditions
% Vcmax_ws maximum carboxylation rate, including the effects of temperature and water stress (mumol/m2leaf/s)
% J_ws     light saturated rate of electron transport, including the effects of temperature and water stress (mumol/m2leaf/s)
% Rd_ws    respiration rate, expressed as a (default) fraction of Vcmax, and thus including the effects of temperature and water stress (mumol/m2leaf/s)
% gamma_star CO2 compensation point, including the effects of temperature and water stress (mol/mol)
% k1_ws    biochemical parameters for the optimization (mumol/m2leaf/s)
% k2_ws    biochemical parameters for the optimization (mumol/m2leaf/s)

%FUNCTIONS CALLED:
%PARAMETERS


%parameters (but see below for the temperature dependences of the metabolic parameters):
PARAMETERS;


%%%%%% CALCULATIONS 

%-------------------- WELL WATERED CASE - Temperature dependences
% %Katul 2000
% Ea=160; %[kJ/mol]
% R=8.314*10^(-3); %[kJ/mol/K]
% alpha_p=0.8;  %unitless
% e_m=0.08;     %mol/mol
% Kc25=300;     %mumol/mol
% Ko25=300;     %millimol/mol
% Coa=210;      %millimol/mol
% Kc=Kc25*exp(0.074*(Tl-25)); %this formulation leads to a slightly lower value than Bernacchi below (and depends in Kc25)
% Ko=Ko25*exp(0.018*(Tl-25));  %this formulation leads to a significantly lower value than Bernacchi below (and depends in Ko25)
% a2_Rub=Kc.*(1+Coa./Ko); %mumol m-2 s-1
% Vcmax_ww=(Vcmax25).*exp(0.088*(Tl-25))./(1+exp(0.29*(Tl-41)));
% gamma=2.6*exp(-0.056*(Tl-25));
% gamma_star=Coa./(2*gamma);
% cp=gamma_star;
% Jmax_ww=Jmax25*exp(Ea*(Tl+273.15-298)./(298*R*(Tl+273.15))); %De Pury and Farquhar 1997


% Vico et al 2013, as extracted from Bernacchi et al 2001
% Rgas = 8.31; %(*universal gas constant (J/(mol K))*)
% Kc = exp(38.05 - 79.43*10^3./(Rgas.*(Tl + 273.15)));%(*% umol/mol (from Bernacchi,2001); the same results, plus or minus 2%, are obtained with eq 5 and 6 in Medlyn 2002 - Part 2*)
% Ko = exp(20.30 - 36.38*10^3./(Rgas.*(Tl + 273.15)));%(*% mmol/mol (from Bernacchi,2001); same as above*)
% Coa = 210;%(*% mmol/mol*)
% gamma_star = exp(19.02 - 37.83*10^3./(Rgas.*(Tl +  273.15)));%(*% ppm (from Bernacchi,2001- these same relationships are used to get the parameters from the A ci curve; the same results are obtained with eq. 12 in Medlyn 2002)*)
% a2_Rub=Kc.*(1+Coa./Ko); %mumol m-2 s-1
% 
% HvV = 116300; Sv = 650; HdV = 202900;
% num = exp(HvV./(Rgas*(298.15)).*(1 - (298.15)./(Tl + 273.15)));
% denomexp = exp((Sv.*(Tl + 273.15) - HdV)./(Rgas.*(Tl + 273.15)));
% Vcmax_ww = Vcmax25*num./(1 + denomexp);
% 
% HvJ = 79500; HdJ = 201000;
% numJ = exp(HvJ./(Rgas*(298.15)).*(1 - (298.15)./(Tl + 273.15)));
% denomexpJ = exp((Sv*(Tl + 273.15) - HdJ)./(Rgas*(Tl + 273.15)));
% Jmax_ww = Jmax25*numJ./(1 + denomexpJ);

%Medlyn 2002 Part 2
%(Medlyn, B. et al: Temperature response of parameters of a biochemically based model of photosynthesis. 
%II. A review of experimental data, Plant Cell Env, 25, 1167-1179, https://doi.org/10.1046/j.1365-3040.2002.00891.x, 2002.
Rgas=8.31;  %(*universal gas constant (J/(mol K))*)
Kc=404.9.*exp(79430.*(Tl-25)./(298.*Rgas.*(Tl+273.15))); %Bernacchi 2001 as reported in Medlyn 2002 (eq. 5)
Ko=278.4.*exp(36380.*(Tl-25)./(298.*Rgas.*(Tl+273.15))); %Bernacchi 2001 as reported in Medlyn 2002 (eq. 6)
gamma_star=42.75.*exp(37830.*(Tl-25)./(298.*Rgas.*(Tl+273.15))); %Bernacchi 2001 as reported in Medlyn 2002 (eq. 6)
Coa = 210;%(*% mmol/mol*)
a2_Rub=Kc.*(1+Coa./Ko); %mumol m-2 s-1

%Vcmax
Ha=70*1000; %J mol^(-1), middle of the observed range (60-80)
Hd=200*1000; %J mol^(-1) - fixed value that gives good results
Topt=41+273.15; %deg K - average value for the two crops
deltaS=Hd/Topt+Rgas*log(Ha/(Hd-Ha)); %inversion of Medlyn 2002 eq 19
numV=1+exp((298*deltaS-Hd)./(298.*Rgas)); 
denomV=1+exp(((Tl+273.15)*deltaS-Hd)./((Tl+273.15).*Rgas)); 
exponV=exp(Ha.*(Tl-25)./(298.*Rgas.*(Tl+273.15)));
Vcmax_ww=Vcmax25.*exponV.*numV./denomV; %peaked function (eq. 17) in Medlyn 2002

%Jmax
Ha=83*1000; %J mol^(-1), average value of the two crops (higher than those for other functional types)
Hd=200*1000; %J mol^(-1) 
Topt=36+273.15; %deg K, average value of the two crops (no clear pattern among species and functional type)
numJ=1+exp((298*deltaS-Hd)./(298.*Rgas)); 
denomJ=1+exp(((Tl+273.15)*deltaS-Hd)./((Tl+273.15).*Rgas)); 
exponJ=exp(Ha.*(Tl-25)./(298.*Rgas.*(Tl+273.15)));
Jmax_ww=Jmax25.*exponJ.*numJ./denomJ; %peaked function (eq. 17) in Medlyn 2002


%------------------- INCLUSION OF WATER STRESS (following Vico and Porporato 2008)
Vcmax_ws=Vcmax_ww*exp(-(-psil/dV)^bV);
phiPSII_ws=phiPSII_ww*exp(-(-psil/dJ)^bJ);
Jphi_ws=1/2*phiPSII_ws*Qp;
J_ws=(Jmax_ww+Jphi_ws-sqrt(Jmax_ww^2+2*Jmax_ww*Jphi_ws+Jphi_ws^2-4*Jmax_ww*Jphi_ws*kphi))/(2*kphi);
Rd_ws=fR*Vcmax_ws;

%-------------------- INPUTS TO OPTIMIZATION MODULE
k1_ws=J_ws/4;
k2_ws=(J_ws/4).*(a2_Rub./Vcmax_ws);

return

