%list of all parameters, and their sources, except those defined locally

%Author: Giulia Vico 
%Date: December 2020
%Part of the codes used to generate the results in the following pubblication. 
%Please refer to that for details on the model, its rationale and parameterization.
%Luan, X. and Vico, G. (2021), Canopy temperature and heat stress are increased by compound high air temperature 
%and water stress, and reduced by irrigation ? A modeling analysis, Hydrol. Earth Syst. Sci., 
%https://doi.org/10.5194/hess-2020-549 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLANT PARAMETERS
%-------- general plant parameters ------
%see Opt_C3_RdgBL_psil for a list of input and default parameters
%leaf features
Vcmax25=83;     %maximum carboxylation rate; mu mol m-2 s-1 - Triticum, based on Vico Porporato 2008/Wulschleger 1993
Jmax25=132;     %maximum electron transport rate; mu mol m-2 s-1 - Triticum, based on Vico Porporato 2008/Wulschleger 1993
fR=0.01;        %fraction of Rd=fR*Vcmax, with Vcmax a function of temperature
RAI=5.6;        %root area index, m2root/m2ground
LAI=2;          %leaf area index, m2leaf/m2ground
dleaf=4/100;    %leaf width (NOT leaf characteristic size), m
hcanopy=0.6;    %canopy height, m
Cd=0.3;         %drag coefficient for the wind profile; Katul et al, Boundary Layer Meteorology, 2004; value for corn and rice
Zr=0.3;         %active rooting depth, m
Tth=30;         %threshold of heat stress damage for a specific phenology, deg C; Triticum, after Saini and Aspinal, 1982


%--------------------- temperature dependences
%see directly metabolic_waterstress.m

%-------------------- response to water stress (after Vico and Porporato 2008)
%Vcmax:
%for Triticum
bV=3.8; %unitless
dV=1.52; %MPa

%Jmax, Jphi and J
%for Triticum
phiPSII_ww=0.59;
bJ=6.2;
dJ=2.6;
kphi=0.9;

%---------------------------marginal water use efficiency and its response to water stress
%the lower the lambda, the higher gs
castar=400; %reference CO2 concentration; mumol/mol
%forbes and grasses, mesic climates, in Manzoni 2011 Functional Ecology (but based on psi leaf)
Lambda_wwstar=981;  %Lambda value under well watered conditions and at ca=400 ppm; mol/mol 
beta0=-1.26;        %1/MPa

%--------------------------cuticular conductance and its dependence on water availability
gcutmax=17.34/10^3;  %residual or cuticolar conductance, under well-watered conditions (based on wheat in Duursma 2019) 
psilstar=-3;         %MPa 
%values of gmin as per Duursma 2019 New Phytol Fig 2c, averages (mmol/m2/s)
%oats: 1.85; peanut: 4.84; cotton: 10.30; sorghum: 11.21; maize: 12.13; soybean: 12.71; millet: 13.08; wheat: 17.34; rice: 19.40


%----------------------------soil to root conductance
%parametrs after Manzoni et al 2013 Advance in Water Resources
dr=0.2/1000; %m root parameter
Edosexponent=8; %exponent of the growth of root with water stress,  if using Edo's root growth (Daly et al 2004, part 1)


%----------------------------root to leaf conductance
gpmax=1.3*10^(-7); %m/(s MPa) %from Vico and Porporato
%for Vico and Porporato 2008
bc=1.1;
dc=1.1; %MPa
%for the alternative approach based on Manzoni
acav=2;
psi50=-2.5; %MPa


%%%%%%%%%%%% FLUIDODYNAMIC PARAMETERS
alphaObukov=0.4/3; %Katul 2004; used only for the first guess of the mixing length - no need to report this
z0=4/1000; %soil roughness length (m), assumed to match that of bare soil (Campbell, Table 5.1), because it 
%is used to impose the (alternative) lower boundary contition on wind speed - NEEDED ONLY IF USING ALTERNATIVE
%BOUNDARY CONDITION



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SOIL PARAMETERS
n=0.43; %soil porosity
%soil features - sandy loam 
psissat=-0.7/10^3; %soil water potential at saturation - WITH SIGN -, MPa
Ksat=0.8/3600/24; %soil hydraulic conductivity at saturation, m/s - enters in the soil to root conductance
bsoil=4.90; %exponent of the soil water retention curve
sfc=0.56; %soil field capacity
%sstar and swilting are not needed, because they 'emerge' from the realtion ET(s) 
s1= 0.62; %was 0.65, but in the paper it is 0.62 - %used only if no diurnal cycle is explored


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ENVIRONMENTAL CONDITIONS (but see also main_plantsolver_optRdgBL_psil_layered.m for those relative to
%precipitation and temperature)
Pair=101; %air pressure, kPa
cairtop=400; %air CO2 concentration, ppm
cloud_frac = 0.1; % cloud fraction (used in layer_conditions.m)
%wind
umeteo=4; %measured wind speed (above the canopy) m/s
zmeteo=2; %height of wind speed measurements (m)


%boundary conditions from the closure scheme
G=0; %start with 0, ie. Tair at the bottom of the canopy is equal to Tsoil; 
%alternative could be to have G=0.1* Rnet (as in Webber 2016, which follows Allen 1998)
%Talow=25;
%Tahi=23;
Ulow=0.01; %m/s - bottom boundary conductance for the wind speed, if imposed as value and not flux 
%it is set to a value greater than 0 to ensure that the air flow is fully turbulent, so that the heat
%exchanges between the air and the ground can be modelled with the turbulent diffisivity, without
%contributions via molecular diffusion (Manzoni 2011 Ecological Modelling - LYCOG)
Uhi=umeteo;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SIMULATION PARAMETERS
layers=500; %number of layers in the whole profile - 400-500 is a bit of a minimum for getting proper profiles
tolerTl=0.1; % toler is acceptable tolerance error for calculating Tl with iterative method (deg C)
tolerpsil=0.001; %toler is acceptable tolerance error for calculating psil with iterative method (MPa);
tolergHa=0.001;  %tolerance for gHa (in aereodynamic conductance) mol/m2/s
maxiter=15; %maximum number of iterations
chopthreshold=1/20; %imaginary to real part ratio below which the imaginary part is neglected
psilguess=[psilstar,0]; %initial guess of the leaf water potential (MPa, with sign); psilstar is the value at which also the cuticular conductance is 0
deltat=0.5; %time step of the model - subdaily scale (hrs): it is immaterial if tmsunrise=tmsunset



%some preliminary calculations
z=[hcanopy];
zM=0.13*hcanopy; %vegetation  momentum roughness length, for dense vegetation (Jones 1992, Ch 3) 
%(used in the aerodynamic_conductance and for the logarithmic wind speed/guess only): 
%must thus match that of the vegetation 3-9 cm according to Campbell Table 5.1; or use the Jones relation
%above
dzero=2/3*hcanopy; %zero plane displacement, for dense vegetation (Jones 1992, Ch 3 has 0.64)







