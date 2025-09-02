function [U,z,uinside,zcanopy,uabove,zabove]=layer_wind_diabatic(zmeteo,umeteo,z,hcanopy,ustar,PsiM)

%Author: Giulia Vico
%Date: September 2020
%This is part of the codes used to generate the results in the following pubblication.
%Please refer to that for details on the model, its rationale and parameterization.
%Luan, X. and Vico, G. (2021), Canopy temperature and heat stress are increased by compound high air temperature
%and water stress, and reduced by irrigation ? A modeling analysis, Hydrol. Earth Syst. Sci.,
%https://doi.org/10.5194/hess-2020-549

%INPUT:
%zmeteo height at which the bulk atmospheric wind speed is defined (i.e., measured) (m)
%umeteo observed/assumed wind velocity at zmeteo (m/s)
%hcanopy canopy height (m)
%z vector containing the layers, till the top of the profile (m)
%ustar friction velocity (m/s)
%PsiM diabatic correction for momemntum transfer (-)

%OUTPUT:
%U:vector containing the wind velocities at each layer in the profiel (m/s)
%z: vector containing the layers, till the top of the profile (m), returned for easy of use
%uinside: wind velocity profile inside the canopy (m/s)
%zcanopy: array of layers, till canopy height (m)
%uabove: wind velocity profile above the canopy (m/s)
%zabove: heights at which the profile above the canopy is calculated (m) - for plotting purposes only

%FUNCTION CALLED:
%PARAMETERS;
%CONSTANTS;

%NOTES
%In dense canopies, Harman and Finnigan (2008) obtained the mean velocity profile as
%U(z)=Uh exp(beta (z-h)/lm)
%with Uh mean velocity at the top of the canopy, beta=ustar/Uh, lm is the mixing length, z has origin at the bottom of the canopy
%Above the canopy, Jones 1992 reports as mean velocity profile
%U(z)=ustar/kvonkarman log((z - dzero)/zM +PhiM) 
%(see also Campbell and Norman 1998 7.24)
%with kvonkarman von kvonkarman constant, ustar the friction velocity, dzero the zero plae displacement and zM the roughness
%height of the vegetation (not of the soil!).
%See Campbell and Norman 1998 Ch 5 for some values or use Jones relation to %canopy height
%Continuity and smoothness in h are imposed; it is further imposed that the velocity at heigh zmeteo is
%umeteo. These conditions allow determining uh, ustar and lm, thus solving the problem

%References: 
%Campbell, G. S., and Norman, J. M.: An introduction to environmental biophysics, Springer, New York City, USA, xv+286 pp., 1998
%Harman, I. N., Finningan, J. J.: Scalar Concentration Profiles in the Canopy and Roughness Sublayer, Boundary-Layer Meteorol, 
%129:323?351, 2008, DOI 10.1007/s10546-008-9328-4
%Jones, H.: Plants and microclimate: a quantitative approach to environmental plant physiology, Cambridge University 
%Press, Cambridge, UK, 428 pp., 1992.
%Siqueira M., Katul G.G.: An Analytical Model for the Distribution of CO2 Sources and Sinks, Fluxes, and Mean Concentration
%Within the Roughness Sub-Layer, Boundary-Layer Meteorol, 2010, DOI 10.1007/s10546-009-9453-8

PARAMETERS;
CONSTANTS;


%--- imposing the conditions above, with no diabatic corrections
%ustar=(kvonkarman*umeteo)/log((-dzero + zmeteo)/zM); %friction velocity (m/s), from imposing the velocity above the canopy
uh=ustar/kvonkarman*(log((hcanopy- dzero)/zM)+PsiM); %wind velocity at the top of the canopy (m/s), from imposing continuity of the profile
beta=ustar/uh; %see Siqueira and Katul 2010
lm=beta*kvonkarman*uh*(-dzero + hcanopy)/ustar; %mixing length, imposing smoothness at z=h

%determine the profile of wind velocity
if length(z)>1
    dz=z(2)-z(1);
    istop=max(find(z<=hcanopy));
    zabove=z(istop+1:end);
    uabove=ustar./kvonkarman*log((zabove - dzero)./zM+PsiM);
    zcanopy=z(1:istop);
    uinside=uh.*exp(beta*(zcanopy - hcanopy)./lm);
    U=uinside;
    U(length(zcanopy)+1:length(zcanopy)+length(zabove))=uabove;    
else %for big leaf models, only the value at the top of the canopy is returned
    U=uh;
    uinside=uh;
    zcanopy=z;
    uabove=uh;
    zabove=z;
end

return