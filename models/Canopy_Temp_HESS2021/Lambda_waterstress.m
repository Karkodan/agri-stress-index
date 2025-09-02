function Lambda_ws=Lambda_waterstress(Lambda_wwstar,ca,psipd)

%Author: Giulia Vico 
%Date: September 2020
%Part of the codes used to generate the results in the following pubblication. 
%Please refer to that for details on the model, its rationale and parameterization.
%Luan, X. and Vico, G. (2021), Canopy temperature and heat stress are increased by compound high air temperature 
%and water stress, and reduced by irrigation ? A modeling analysis, Hydrol. Earth Syst. Sci., 
%https://doi.org/10.5194/hess-2020-549 

%INPUT
% Lambda_wwstar (mol/mol): marginal water use efficiency under well-watered conditions
% ca (mumol/mol): atmospheric CO2 concentration
% psipd (Mpa): predawn soil water potential (with sign)

%OUTPUT
%Lambda_ws: marginal water use efficiency, including the effects of water stress (mol/mol)

%FUNCTIONS CALLED:
%PARAMETERS

%NOTES: 
% This follows Manzoni et al 2011 Functional Ecology and Zhou 2013 Agr and Forest Met BUT using values recalculated by 
% Manzoni (Aug 2018 personal communication), so as to consider the dependence on psi on predawn leaf water potential 
% as opposed to psileaf directly

%Note that, according to Manzoni et al, if the cuticular conductance is neglected
%a parabolic dependence should be used; conversely, Zhou simply uses an exponential dependence. 
%However, the two dependences are basically coincident at moderate water stress. 
%This is the most common condition for crops and, as such, the choice of the dependence of lambda on psi leaf
%bears small consequences for our applications. The exponential dependence is preferred because it
%has just one parameter and, as such, its fitting to the data is more robust.
%Between the two formulations (Manzoni's and Zhou's), that of Manzoni was chosen. 

%References
%Manzoni, S., Vico, G., Katul, G., Fay, P. A., Polley, W., Palmroth, S., and Porporato, A.: Optimizing stomatal 
%conductance for maximum carbon gain under water stress: a meta?analysis across plant functional types and climates, 
%Funct Ecol, 25, 456-467, https://doi.org/10.1111/j.1365-2435.2010.01822.x, 2011.
%Zhou, S., Duursma, R. A., Medlyn, B. E., Kelly, J. W., and Prentice, I. C.: How should we model plant responses 
%to drought? An analysis of stomatal and non-stomatal responses to water stress, Agr Forest Met, 182, 204-214,
%https://doi.org/10.1016/j.agrformet.2013.05.009, 2013.


PARAMETERS;

%calculations
Lambda_ws=Lambda_wwstar*ca/castar*exp(beta0*psipd);


return