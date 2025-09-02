%Authors: Xiangyu Luan, Giulia Vico 
%Date: December 2020
%Part of the codes used to generate the results in the following pubblication. 
%Please refer to that for details on the model, its rationale and parameterization.
%Luan, X. and Vico, G. (2021), Canopy temperature and heat stress are increased by compound high air temperature 
%and water stress, and reduced by irrigation ? A modeling analysis, Hydrol. Earth Syst. Sci., 
%https://doi.org/10.5194/hess-2020-549 

%KEY OUTPUTS:
% The run returns the time series of the key model outputs, including
% s soil moisture (unitless)
% Tl canopy temperature (deg C)
% Anleaf leaf-level photosynthetic rate (mu mol/m^2/s)
% ciopt CO2 concetration at the photosynthetic site, under optimal stomatal conductance (mol/mol)
% Lambda_ws marginal water use efficiency (mol/mol)
% Vcmax maximum carboxylation rate, considering canopy temperature and water status (mu mol/m^2leaf/s)
% ETcanopy evapotranspiration rate (per unit ground area) (mol/s/m2ground)
% gal, gbl,gh, gha, gp, gs, gbl gscopt, gsrp, gv conductances (mol/m^2/s)
% psil psir psis water potentials of leaf, root and soil respectively (MPa)
% VPD canopy water pressure decificit (mol/mol)
% rainsubd daily precipitation (mm)
% T_airsubd daily maximum air temperature (deg C)
% Rabs absorbed radiation (W/m^2)
% Q_Short and Q_Long incomning short and longwave radiation (W/m^2)

%FUNCTIONS CALLED:
%CONSTANTS
%PARAMETERS
%weathergenerator
%aereodynamic_cond
%layer_wind_diabatic
%Opti_C3_RdgBL_psil, which calls:
    %metabolic_waterstress 
    %Lambda_waterstress 
    %gsAciopt_RdgBL 
    %soil2leaf_conductance 

%metabolic_waterstress
%Tleaf_layers
%bigleaf_long_shortwave, which calls
    %solar_dec_azimuth
%soilbalanceET  



close all;
clc;
clear;

%-------- KEY MODEL SETTINGS
%-- initial condition for soil moisture
s_int=0.56; %initial soil moisture (at the beginning of the simulation)
initialcond=0; %if 1, then soil moisture starts at the same value, s_int, at the beginning of each replicate;
%otherwise, the last soil moisture of the previous replicate is used, so that, to all effects, the initial
%conditions of each replicate (except the first one) is stochastic
%-- irrigation
irrigation=1; %if 1, irrigation is employed
diagnosticopt=0; %defines whether the optimization module produce a plot showing the convergence; if negative, no plot is set
%-- number of days simulated
dds=21; %duration of the simulation (days)
repls=2; %number of replicated simulations (501 in Luan and Vico 2021 HESS)


%-------- DEFINITION OF CONSTANTS AND PARAMETERS
CONSTANTS;
PARAMETERS;


%-------- DEFINITION OF ENVIRONMENTAL CONDITIONS (not set in PARAMETERS)
%-- air humidity
RH=40; %relative humidity (in percent)
%-- precipitation parameters
rainY=1000; %annual rainfall (mm)
Fr=0.45;    %rainfall frequency (1/d)
%-- air temperature
Tmaxmean=24.0; %average of air temperature distribution, unit degree Celsius - to be interpreted as daily maximum
k1=0.81;  %drift parameter of air temperature (d).
k3=32.8;  %diffusion parameter of air temperature.
%-- irrigation
if irrigation==1
    stilde=0.25; % soil moisture corresponding to the intervention point for irrigation
else
    stilde=0;
end
shat=sfc+0.1;   % soil moisture corresponding to the target level of irrigation
%-- solar radiation
R_incnoon=1000; %W m^(-2)
%-- location
firstDOY=210; %first day of the year (doy) of the simulation
LAT=45; %latitude of the location (in decimal degrees)
LONG=10; %ongitude of the location (in decimal degrees) - it does not bear any consequence
tm=12; %time of the day - noon only (one time point per day)


%----------- PRESET VECTORS
timestepsday=1;
ts=[1:1:dds*repls]; %days
%parameters for the daily integration of ET losses
%ET must be determined from gv at noon (cumulated over all layers), assuming a parabolic dependence on time of the day
duration=12; %hrs per day during which photosynthesis and transpiration occur
coeff=2/3; %assuming the calculated ET is the max daily value, and a parabolic temporal evolution over teh day; the length of the day is immaterial for this constant
%variables needing one value per day
[psipd,sday]=deal(zeros(1,dds*repls));
%variables needing one value per time step
[psis,psir,s,ETcanopy,gal]=deal(zeros(1,dds*timestepsday*repls));
%variables needed one value per time step and layer inside the canopy only
[Anleaf,Eleaf,Eplant,ciopt,gscopt,gsrp,gsr,gp,gs,gsbl,gbl,gv,gh,gha,Q_Short,Qp,Q_Long,Tl,Tlnew,errTl,psil,Vcmax_ws,Lambda_ws]=deal(zeros(length(z),dds*timestepsday*repls));
%variables needed one value per time step and layer over the whole profile
[U,Uexpln,T_air,es,VPD,cair]=deal(zeros(length(z),dds*timestepsday*repls));
%initialization of soil moisture (at least for the first replicate)
s(1)=s_int;sday(1)=s(1);


%--------------------- GENERATION OF THE RANDOM WEATHER CONDITIONS
%%% these lines are to be substituted with data retrieval commands, if observations are used
%weathergenerator generates daily rainfall and daily max temperature - one value per day
dt=1; %days
Tin=Tmaxmean; %initial air temperature value for the first day simulation
[rain,Tmeteo]=weathergenerator(dt,Tmaxmean,Tin,rainY,Fr,dds*repls,k1,k3);


%--------------------- SPECIAL CONDITIONS FOR TESTING - this overrides the above generated conditions
% comment/uncomment specific blocks as needed
%-- 1) dry down
%rain=rain*0;
%Tmeteo=Tmeteo*0+Tmaxmean;
%-- 2) only air temperature changing (but well watered conditions)
% rain=rain*0+10;
% Tmeteo=[15:1:15+dds-1]; %deg C
%-- 3) parameters to compare with Campbell and Norman - Example 14.1
% R_incnoon=1100; %W/(m2 s)
% Tmeteo=Tmeteo*0+38; %deg C
% dleaf=10/100; %m
% umeteo=3.5; %m/s


%%%%%%%%%%%%%%%%
%CALCULATIONS
%%%%%%%%%%%%%%%%


%FURTHER PREPARATION OF ENVIRONMENTAL CONDITIONS AND PLOTTING
Uhi=umeteo;
%solar radiation and air temperature at measurement height
R_inc=R_incnoon+zeros(1,dds*repls);
T_airsubd=Tmeteo+zeros(length(z),1);
cair(1:length(z),1:dds*repls)=cairtop+zeros(length(z),dds*repls);
rainsubd=rain;

%VPD at measurement height
estopsubd=a.*exp((b.*T_airsubd(end,:))./(T_airsubd(end,:)+c)); %air vapor pressure at saturation [kPa, with temp in deg C]
Dtopsubd=estopsubd*(1-RH/100); %vapor pressure deficit, kPa
VPDtopsubd=Dtopsubd/Pair; %converts D (KPa) to  VPD  (mol water/mol air)

%-- plotting of environmental conditions used for the simulations
figure(4)
subplot(4,1,1)
plot(ts, R_inc)
xlabel('time');ylabel('R_{inc} (W m^{-2})')
subplot(4,1,2)
plot(ts,T_airsubd)
xlabel('time');ylabel('T_{air}(z=z_{m}) (^{o}C)')
subplot(4,1,3)
plot(ts,VPDtopsubd)
xlabel('time');ylabel('D (z=z_{m}) (mol mol^{-1})')
subplot(4,1,4)
plot(ts,rainsubd)
xlabel('time');ylabel('Precipitation (mm)')


%DETERMINATION OF CANOPY TEMPERATURE

i=1;d=1;id=1;
while d<=dds*repls %cycle on day - NOTE: if gs becomes extremely low, the loop is exited; this is to avoid running into numerical problems in the optimization cycle (where the demand for water becomes negative)
    %initialization of soil moisture
    if floor((d-1)/dds)==(d-1)/dds && initialcond==1 %i.e., if the initial soil moisture for each replicate is set to be the same and equal to s_int; otheerwise the end of the previous run is sused
        s(id)=s_int;sday(d)=s_int;
    end
    waitbar(d/dds*repls)
    psipd(d)=psissat*sday(d)^(-bsoil); %predawn soil water potential (MPa) - explicitly distinguished from instantaneous psisoil for easy extension to sub-daily scale (lambda does not respond instantaneously to leaf water potential)
    DOY=firstDOY-1+d;
    psis(i)=psissat*s(i)^(-bsoil); % soil water potential (MPa)
    Tl(:,i)= T_airsubd(:,i); %APPROXIMATION
    %bigleaf_long_shortwave provides short, PAR and longwave radiation;
    %NOTE: longwave radiation is actually not used other than to calculate Rabs for the plotting and it is based on whatever is the
    %current estimate of Tleaf (no binomial expansion + linearization used)
    [Q_Short(:,i),Qp(:,i),Q_Long(:,i)]=bigleaf_long_shortwave(LAT, R_inc(:,i), T_airsubd(:,i), Tl(:,i), DOY, tm);
    if max(Q_Short(:,i))<10 %handling of very low (unrealistically?) values of solar radiation, to avoid numerical issues
        Rabs(:,i)=Q_Short(:,i)+Q_Long(:,i); %total absorbed radiation
        psil(:,i)=zeros(1:length(z),1)+psis(i);
        psir(:,i)=zeros(1:length(z),1)+psis(i);
        gs(1:length(z),i)=zeros(1:length(z),1)+max(-gcutmax/psilstar*psil(:,i)+gcutmax,0);
        U(:,i)=Uhi; %APPROXIMATION
        [gHa,HH,ustar,PsiM,PsiH]=aereodynamic_cond(T_airsubd(end,i),Tl(1,i));
        gal(:,i)=gHa/LAI;
        gbl(:,i)=1.4*0.147*sqrt(U(:,i)./(0.7*dleaf));      %leaf boundary layer conductance to water vapor
        gh(:,i)=1.4*0.135*sqrt(U(:,i)./(0.7*dleaf));       %mol/m2/s  %leaf boundary layer conductance to heat
        gsbl(i)=gs(i)*gbl(i)/(gs(i)+gbl(i));               %per unit leaf area
        gv(i)=gsbl(i)*gal(i)/(gsbl(i)+gal(i));             %per unit leaf area
        es(:,i) = a.*exp((b.*T_airsubd(:,i))./(T_air(:,i)+c)); %air vapor pressure at saturation (kPa, with temp in deg C)
        VPD(:,i)=es(:,i)*(1-RH./100)/Pair;                 %vapor pressure deficit (mol H2O/mol air)
        [~,~,Rd_ws,~,~,~]=metabolic_waterstress(Vcmax25, Jmax25,psil(i),0,Tl(i)); %only respiration matters, at such low solar radiation
        Anleaf(i)=-Rd_ws;
        deltae=VPD(:,i); %deltae=VPD because Tl assumed to be equal to Tair under so low radiation
        gsla=gsbl(i)*gal/(gsbl(i)+gal);
        Eleaf(i)=gsla*deltae;
        Eplant(i)=Eleaf(i)*18*10^(-3)*10^(-3);
    else %when radiation conditions are similar to those typical of the most of the day
        psis(i)=psissat*s(i)^(-bsoil); % soil water potential (MPa)
        
        %-- definition of the initial guesses for the scalar and momentum transport equations
        Uprevious=Uexpln(:,i)'; Taprevious=-1;          %set to -1, because a constant Tair will stop the simulation
        gbl(:,i)=1.4*0.147*sqrt(U(:,i)./(0.7*dleaf));   %leaf boundary layer conductance to water vapor
        gh(:,i)=1.4*0.135*sqrt(U(:,i)./(0.7*dleaf));    %leaf boundary layer conductance to heat (mol/m2/s)
        
        %-- DETERMINATION OF CANOPY TEMPERATURE (iteratively)
        iternumber=0; errTl(1,i)=2*tolerTl; %high errTl to enter the first round of the loop
        while max(errTl(:,i))>tolerTl && iternumber<maxiter %errTl(i) is the difference in Tl between two subsequent iterations; if small, it exists the loop
            %air temperature at measurement height
            T_air(:,i)=T_airsubd(end,i);
            %-- recalculate wind profile (and particularly wind at the top of the canopy), including diabatic corrections
            [gHa,HH,ustar,PsiM,PsiH]=aereodynamic_cond(T_air(end,i),Tl(1,i));
            [Uexpln(:,i),~,~,~,~,~]=layer_wind_diabatic(zmeteo,umeteo,z,hcanopy,ustar,PsiM);
            U(:,i)=Uexpln(:,i);
            gal(:,i)=gHa/LAI; %the aereodynamic conductance in the bulk atmosphere is the same for all scalars;
            %it is divided by LAI to go from mol/m2 ground/s to mol/m2 leaf/s, to match the fact that it is
            %later put in series with the stomatal conductance, not canopy condcuctance
            
            %based on the calculated profiles of U and Ta, define or redefine some additional profiles
            es(:,i) = a.*exp((b.*T_air(:,i))./(T_air(:,i)+c));  %air vapor pressure at saturation (kPa, with temp in deg C)
            VPD(:,i)=es(:,i)*(1-RH./100)/Pair;                  %vapor pressure deficit (mol H2O/mol air)
            gbl(:,i)=1.4*0.147*sqrt(U(:,i)./(0.7*dleaf));       %leaf boundary layer conductance to water vapor
            gh(:,i)=1.4*0.135*sqrt(U(:,i)./(0.7*dleaf));        %leaf boundary layer conductance to heat (mol/m2/s)
            gha(:,i)=gal(:,i)*gh(:,i)./(gal(:,i)+gh(:,i));      %both expressed on a per unit leaf area basis (see below for the upscaling)
            
            %-- determination of canopy activity, based on metabolic functioning and stomatal limitations
            [Anleaf(i),Eleaf(i),Eplant(i),gscopt(i),ciopt(i),psil(i), ...
                Vcmax_ws(i),Lambda_ws(i),gsrp(i),gsr(i),gp(i)]= ...
                Opti_C3_RdgBL_psil...
                (Vcmax25,Jmax25,Lambda_wwstar,Tl(i),psilguess,psipd(d),psis(i),Qp(i),cair(i),...
                VPD(i),T_air(i),gbl(i),gal(i),psissat,Ksat,bsoil,Zr,RAI,LAI,diagnosticopt);
            if abs(imag(gscopt(i)))/abs(real(gscopt(i)))<chopthreshold
                gs(i)=real(gscopt(i));
            else
                gs(i)=gscopt(i);
            end
            %%Uncomment the following lines to get the step by step results, for debugging
            %                     windspeed=U(i)
            %                     optimalstcond=gscopt(i)
            %                     leaftemperature=Tl(layer, i)
            %                     soilpredawnwaterpot=psipd(d)
            %                     soilwaterpot=psis(i)
            %                     PAR=Qp(i)
            %                     airCO2conc=cair(i)
            %                     vaporpressuredeficit=VPD(i)
            %                     maxcarboxylationrate=Vcmax_ws(i)
            %                     marginalWUE_waterstress=Lambda_ws(i)
            %                     psileaf= psil(i)
            %                     cio=ciopt(i)
            %                     El=Eleaf(i)
            %                     pause
            
            %-- calculation of canopy temperatute: here Tair is interpreted as the bulk atmospheric temperature,
            %so i) also the aereodynamic conductance is needed; ii) everything can be expressed on a per unit ground area
            gsbl(i)=gs(i)*gbl(i)/(gs(i)+gbl(i));        %per unit leaf area
            gv(i)=gsbl(i)*gal(i)/(gsbl(i)+gal(i));      %per unit leaf area
            [Tlnew(i)]=Tleaf_layers(Pair,LAI*gha(i),T_air(i),Pair*VPD(i),LAI*gv(i),es(i),Q_Short(i),LAI); 
            %In the input to Tleaf_layers, the LAI in front of gv transforms the series gs, gbl, gal from a per unit leaf area to a per ground area basis: as such this is a canopy conductance
            errTl(i)=abs(Tlnew(i)-Tl(i));
            Tl(i)=Tlnew(i);
            iternumber=iternumber+1;
        end %end of while cycle, calculating iteratively all the quantities
        
        %recalculate the radiations, based on the final leaf temperature (for plotting purposes)
        [Q_Short(:,i),Qp(:,i),Q_Long(:,i)]=bigleaf_long_shortwave(LAT, R_inc(:,i), T_airsubd(:,i), Tl(:,i), DOY, tm);
        Rabs(:,i)=Q_Short(:,i)+Q_Long(:,i); %total absorbed radiation
        psir(i)=(gsr(i)*psis(i)-gsrp(i)*(psis(i)-psil(i)))/gsr(i);
        
    end %end of the day vs night cycle / low light vs normal light cycle
    
    %upscaling of water fluxes at the whole canopy level and update of soil moisture balance
    %calculation of the losses via ET
    
    %slope of the saturation vapour pressure-temperature relationship, calculated at the intermediate temp between Ta and Tl
    Tinterm=(Tl(1:length(z),i)+T_air(length(z),i))/2;
    ss=4098.*(0.6108.*exp(17.27.*Tinterm./(Tinterm+237.3)))./(Tinterm+237.3).^2; %kPa/degC (NB This formula is slightly different from that reported in Jones, Appendix 4, but the numerical results are by and large the same, when considering this is in kPa/deg C
    
    %humidity gradient, driving the transpiration
    deltae=VPD(1:length(z),i)+ss./Pair.*(Tl(1:length(z),i)-T_air(1:length(z),i)); %Pair is the air pressure in kPa, so that deltae is in mol/mol, NOT in Pa (or kPa)
    
    %assuming constant leaf area density, it is enough to average the contributions of the different layers
    isok=isfinite(deltae); %should there be any NaN in Tl, then they are not considered (if considered ETcanopy=NaN, which then leads to a ds=NaN, which turns in snew=smax)
    ETcanopy(i)=mean(deltae(isok).*gv(isok,i))*LAI; % ET is in mol/s/m2ground (gv was on a per unit leaf area basis)
    [s(i+1),ds]=soilbalanceET(s(i),ETcanopy(i),rainsubd(i),duration,coeff,sfc,s1,Zr,n,Ksat,bsoil,stilde,shat);
    i=i+1;id=id+1;
    sday(d+1)=s(i); %updating of the daily soil moisture (to be used for the predawn soil water potential)
    d=d+1;
end

Max=Tl;
for i=1:dds
    if Max(i)==0
        disp('There is zero value in Max')
    end
end
% MeanMax=mean(Max); % Return max leaf temperature of all layers for the whole simulation time series.

HTls=find(Max>Tth); % Find High Temperature of leaf which exceeds plant threshold value

clear a a_L a_NIR a_PAR a_S aa acav alphaObukov 
clear b bb bc beta0 bJ bsoil bV c
clear cairtop castar CF chopthreshold cloud_frac col1 Cp Cd
clear d dc deltae deltat diagnosticopt dJ dleaf DOY dr ds dt duration dV dzero
clear Edosexponent epsil epssoil  Eplant
clear fR G gcutmax gg gHa 
clear h h1 HH HTls id initialcond irrigation isol iternumber 
clear Jmax25 k_C k_H k_W kd kphi Ksat kvonkarman
clear lambda_wwstar lambdKg layer1 layers led maxiter n NIR_frac P
clear PAR_frac phiPSII_ww PM_air PM_w psi50 PsiH psilguess psilstar PsiM psissat q psipd
clear rain R_incnoon Rgas rho_d_NIR rho_d_PAR rho_molgas rho_w s_int sday sigm sign_NIR sigm_PAR ss
clear Taprevious temp timestepsday Tin Tinterm Tlnew tm Tmeteo tolergHa tolerpsil 
clear Uhi Ulow Uprevious ustar Vcmax25 z0 zM zmeteo

%%%%%%%%%%% FIGURES
%-------------------- figures for quick checks - plots the whole time series (including replicates)
temp=find(s>0);
isok=temp(1:end-1);%s=0 means that the numerical scheme has not worked properly (typically at low soil moisture availability,
%the optimization scheme returns at some point a negative gs; to avoid error messages, any gs<5*10^(-3) will
%force the cycle to exit

layer1=1;
col1=[0.2,0.8,0.2];

figure(1)
subplot(5,4,1)
plot(ts,rainsubd(1:isok(end)),'-','Color',col1)
ylabel('Precipitation (mm)');xlabel('t (d)')

subplot(5,4,5)
plot(ts,T_airsubd(1,1:isok(end)),'-','Color',col1);
ylabel('T_{air} (^oC)'); xlabel('t (d)')

subplot(5,4,9)
plot(ts,VPD(layer1,1:isok(end)),'-','Color',col1);
ylabel('D_l (mol mol^{-1})');xlabel('t (d)')

subplot(5,4,13)
plot(ts,Rabs(layer1,1:isok(end)),'-','Color',col1);hold on;
plot(ts,Q_Short(layer1,1:isok(end)),':','Color',col1);hold on;
plot(ts,Q_Long(layer1,1:isok(end)),'--','Color',col1);hold off;
ylabel('Radiation (W m^{-2})');xlabel('t (d)')
legend('R_{abs}','Q_{short}','Q_{long}');

subplot(5,4,2)
plot(ts,psil(layer1,1:isok(end)),'-','Color',col1); hold on
plot(ts,  psir(1:isok(end)),':','Color',col1);
plot(ts,psis(1:isok(end)),'-','Color',[0.5,0,0.0]);hold off
legend('\psi_{l}','\psi_r','\psi_{s}')
ylabel('\psi (MPa)');xlabel('t (d)')

subplot(5,4,6)
plot(ts,s(1:isok(end)),'-','Color',[0.5,0,0.0]);hold off
ylabel('s (-)');xlabel('t (d)')

subplot(5,4,10)
plot(ts,ETcanopy(1:isok(end))*18/10^6*duration*3600*coeff*10^3,'-','Color',[0.5,0,0.0]);hold off
ylabel('ET (mm d^{-1})');xlabel('t (d)')

subplot(5,4,14)
plot(s(1:isok(end)),100*ETcanopy(1:isok(end))/max(ETcanopy(1:isok(end))),'-','Color',col1);hold on
plot(s(1:isok(end)),100*Anleaf(layer1,1:isok(end))/max(Anleaf(layer1,1:isok(end))),':','Color',col1);hold off
legend('ET','A_{net}')
ylabel('% of max (-)');xlabel('s (-)')

subplot(5,4,18)
plot(gs(layer1,1:isok(end)),ETcanopy(1:isok(end))/max(ETcanopy(1:isok(end))),'-','Color',col1);hold on
plot(gs(layer1,1:isok(end)),Anleaf(layer1,1:isok(end))/max(Anleaf(layer1,1:isok(end))),':','Color',col1);hold off
legend('ET','A_{net}')
ylabel('% of max (-)');xlabel('g_{sc} (mol m^{-2} s^{-1})')

subplot(5,4,3)
plot(ts,gs(layer1,1:isok(end)),'-','Color',col1)
ylabel('g_{sc} (mol m^{-2} s^{-1})');xlabel('t (d)')

subplot(5,4,7)
plot(ts,Anleaf(layer1,1:isok(end)),'-','Color',col1)
ylabel('A_{net} (\mu mol m^{-2} s^{-1})'); xlabel('t (d)')

subplot(5,4,11)
plot(ts,gsrp(layer1,1:isok(end)),'-','Color',[0.5,0,0.0]);hold on
plot(ts,gsr(layer1,1:isok(end)),':','Color',[0.5,0,0.0]);hold on
plot(ts,gp(layer1,1:isok(end)),'--','Color',[0.5,0,0.0]);hold on
hold off
legend('gsrp','gsr','gp'); ylim([0 1.5*gpmax])
ylabel('g (m s^{-1} MPa{-1})'); xlabel('t (d)')

subplot(5,4,4)
plot(ts,Tl(layer1,1:isok(end)),'-','Color',col1)
ylabel('T_{l} (MPa)');xlabel('t (d)')

subplot(5,4,8)
plot(T_airsubd(1:isok(end)),Tl(layer1,1:isok(end)),'.b','MarkerFaceColor',col1);hold on
plot([min(T_airsubd)*0.9,max(T_airsubd)*1.1],[min(T_airsubd)*0.9,max(T_airsubd)*1.1],':k'); hold off
xlabel('T_a (^oC)');ylabel('T_c (^oC)')

subplot(5,4,12)
plot(T_airsubd(1:isok(end)),Tl(layer1,1:isok(end))-T_airsubd(1:isok(end)),'.b','MarkerFaceColor',col1);hold on;
plot([min(T_airsubd)*0.9,max(T_airsubd)*1.1],[0,0],':k'); hold off
xlabel('T_a (^oC)');ylabel('T_c-T_a (^oC)')

subplot(5,4,16)
plot(s(1:isok(end)),Tl(layer1,1:isok(end))-T_airsubd(1:isok(end)),'.b','MarkerFaceColor',col1);hold on;
plot([min(s(1:isok(end)))*0.9,max(s(1:isok(end)))*1.1],[0,0],':k'); hold off
xlabel('s (-)');ylabel('T_c-T_a (^oC)')

subplot(5,4,15)
plot(ts,real(Vcmax_ws(layer1,1:isok(end))),'-','Color',col1)
ylabel('V_{c,max} (\mu mol m^{-2} s^{-1})'); xlabel('t (d)')

subplot(5,4,19)
plot(ts,real(Lambda_ws(layer1,1:isok(end))),'-','Color',col1)
ylabel('\lambda_{w} (mol mol^{-1})'); xlabel('t (d)')

figure(2)
subplot(2,2,1)
plot(gs(layer1,1:isok(end)),Tl(layer1,1:isok(end))-T_airsubd(1:isok(end)),'.','MarkerFaceColor',col1)
xlabel('g_s (mol m^{-2} s^{-1})');ylabel('T_l-T_a(^oC)');
subplot(2,2,2);
plot(s(1:isok(end)),Tl(layer1,1:isok(end))-T_airsubd(1:isok(end)),'.','MarkerFaceColor',col1);
h=refline(0,0);
ylim([-10,10]);
set(h,'Color',[0.5,0,0.0],'LineWidth',1);
xlabel('Soil moisture (-)');ylabel('T_l-T_a (^oC)');
subplot(2,2,3);
plot(T_airsubd(1:isok(end)),Tl(layer1,1:isok(end)),'.','MarkerFaceColor',col1);
h=refline(1,0);
h1=refline(0,34);
set(h,'Color','b','LineWidth',1);
set(h1,'Color','b','LineWidth',1);
xlabel('T_a (^oC)');ylabel('T_c (^oC)');
subplot(2,2,4)
fig = figure(2);
left_color = [.5 0 0];
right_color = [0 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
plot(s(1:isok(end)),'Color',[0.5,0,0.0])
xlabel('time (day)')
ylabel('Soil moisture s (-)')
yyaxis right
plot(T_airsubd,':r'); hold on
plot(Tl(layer1,1:isok(end)),'Color',col1);
ylabel('T (^oC)')
legend('s','T_a','T_c')



fig = figure(3);
left_color = [.5 0 0];
right_color = [0 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
plot(s(1:isok(end)),'-.','Color',[0.5,0,0.0],'LineWidth',7)
xlabel('Time (day)','Fontsize',40)
ylabel('Soil Moisture s (-)','Fontsize',40)
yyaxis right
plot(T_airsubd,':r','LineWidth',7); hold on
plot(Max,'color',col1,'LineWidth',7);
hold off
ylabel('T (^oC)','Fontsize',40)
led=legend('s','T_{a}','T_{c}');
set(gca,'fontsize',40,'looseInset',[0 0 0 0]);
set(led,'fontsize',40);





figure(1)
