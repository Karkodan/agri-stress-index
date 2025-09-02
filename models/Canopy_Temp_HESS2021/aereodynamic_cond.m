function [gHa,H,u_star,PsiM,PsiH]=aereodynamic_cond(Tzmeteo,Tc)

%Author: Giulia Vico
%Date: September 2020
%This is part of the codes used to generate the results in the following pubblication.
%Please refer to that for details on the model, its rationale and parameterization.
%Luan, X. and Vico, G. (2021), Canopy temperature and heat stress are increased by compound high air temperature
%and water stress, and reduced by irrigation ? A modeling analysis, Hydrol. Earth Syst. Sci.,
%https://doi.org/10.5194/hess-2020-549


%This function determines the conductance between the 'canopy surface' and a generic height zmeteo above the canopy, 
%where measurements are taken.
%The code follows Campbell and Norman 1998, Ch 7, sections 7.5 and 7.6.
%The same approach is used in Webber et al 2016 Environmental modelling and software. 
%Instead Tuzet et al 2003 Plant Cell and Environment follows a similar approach, but considers an additional conductance,
%going from hs (which appears to be the 'functional height' at which the big leaf is located - unclear how defined) and hc (which
%appears to be the actual canopy height. But then the conductances are summed up, so the end result will be
%the same.
%Some iterations are performed directly within the code below.

%Additional references
%Bonan, G.: Climate change and terrestrial ecosystem modeling, Cambridge University Press, Cambridge, UK, xx+438 pp., 2019
%Campbell, G. S., and Norman, J. M.: An introduction to environmental biophysics, Springer, New York City, USA, xv+286 pp., 1998
%Katul, G.G., Konings, A.G., Porporato A.: Mean Velocity Profile in a Sheared and Thermally Stratified
%Atmospheric Boundary Layer, Physical Review Letters, 107, 268502, 2011
%Tuzet, A., Perrier, A., and Leuning, R.: A coupled model of stomatal conductance, photosynthesis and transpiration, 
%Plant Cell Env, 26, 1097-1116, https://doi.org/10.1046/j.1365-3040.2003.01035.x, 2003

%INPUT
%Tzmeteo: air temperature at height zmeteometeo (i.e., measured above the canopy), deg C
%Tc: canopy temperature, deg C

%OUTPUT
%gHa: aereodynamic conductance (mol/m2 ground/s)
%H: sensible heat fux (J/m2/s)
%u_star: friction velocity (m/s)
%PsiM: diabatic correction factors for momentum transfer
%PsiH: diabatic correction factors for heat transfer


%FUNCTIONS CALLED
CONSTANTS;
PARAMETERS;


zH=0.2*zM;    %roughness length for heat (and other scalars): follows Campbell, eq. 7.19; 
              %zM is the momentum roughness length for the vegetation, and it is one parameter

%checks and adjusts units of temperature
if Tzmeteo<270
        TaK=Tzmeteo+273.15;
else
        TaK=Tzmeteo;
        Tzmeteo=Tzmeteo-273.15;
end
if Tc>270
    Tc=Tc-273.15;
end

%initial guesses just to enter in the cycle the first time
PsiMguess=-0.9;
PsiM=PsiMguess;
PsiH=0.6*PsiMguess;
gHa=kvonkarman^2*rho_molgas*umeteo/((log((zmeteo-dzero)/zM)+PsiM)*(log((zmeteo-dzero)/zH)+PsiH));
gHanew=gHa+2*tolergHa; 
iter=0;
while abs(gHanew-gHa)>tolergHa && iter<=maxiter
    H=kvonkarman^2*rho_molgas*Cp*umeteo*(Tc-Tzmeteo)/((log((zmeteo-dzero)/zM)+PsiM)*(log((zmeteo-dzero)/zH)+PsiH));
    u_star=umeteo*kvonkarman/(log((zmeteo-dzero)/zM)+PsiM);
    zmeteoeta=-kvonkarman*gg*(zmeteo-dzero)*H/(rho_molgas*Cp*TaK*u_star^3); 
    %Campbell and Norman 1998 has just zmeteo; Katul et al 2011 PRL the same, but in the absence of vegetation; 
    %most of the references including vegetation (e.g. Bonan 2019 Eq. 6.30) have zmeteo-dzero
    %DIABATIC CORRECTION FACTORS (conditions are adiabatic only when sensible heat flux H=0)
    if H>0
        %assuming unstable flow, i.e. normal diurnal conditions, with the canopy warmer than the air above it
        PsiH=-2*log((1+(1-16*zmeteoeta)^(1/2))/2);
        psiM=0.6*PsiH;
    else
        %for stable flow
        PsiH=6*log(1+zmeteoeta);
        PsiM=PsiH;
    end
    gHanew=kvonkarman^2*rho_molgas*umeteo/((log((zmeteo-dzero)/zM)+PsiM)*(log((zmeteo-dzero)/zH)+PsiH));
    iter=iter+1;
end

gHa=gHanew;
H=kvonkarman^2*rho_molgas*Cp*umeteo*(Tc-Tzmeteo)/((log((zmeteo-dzero)/zM)+PsiM)*(log((zmeteo-dzero)/zH)+PsiH));
u_star=umeteo*kvonkarman/(log((zmeteo-dzero)/zM)+PsiM);

return



