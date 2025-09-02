function [Anleaf, Eleaf, Eplant, gscopt, ciopt, psil,Vcmax_ws,Lambda_ws,gsrp,gsr,gp] = Opti_C3_RdgBL_psil...
    (Vcmax25,Jmax25,Lambda_wwstar,Tl,psilguess,psipd,psis,Qp,ca,VPD,Ta,gbl,ga,psissat,Ksat,bsoil,Zr,RAI,LAI,diagnostic)

%Author: Giulia Vico
%Date: September 2020
%This is part of the codes used to generate the results in the following pubblication.
%Please refer to that for details on the model, its rationale and parameterization.
%Luan, X. and Vico, G. (2021), Canopy temperature and heat stress are increased by compound high air temperature
%and water stress, and reduced by irrigation ? A modeling analysis, Hydrol. Earth Syst. Sci.,
%https://doi.org/10.5194/hess-2020-549

%Key additional references:
%Manzoni, S., Vico, G., Katul, G., Fay, P. A., Polley, W., Palmroth, S., and Porporato, A.: Optimizing stomatal
%conductance for maximum carbon gain under water stress: a meta?analysis across plant functional types and climates,
%Funct Ecol, 25, 456-467, https://doi.org/10.1111/j.1365-2435.2010.01822.x, 2011.

%Vico, G., Manzoni, S., Palmroth, S., Weih, M., and Katul, G.: A perspective on optimal leaf stomatal conductance
%under CO2 and light co-limitations, Agr Forest Met, 182, 191-199, https://doi.org/10.1016/j.agrformet.2013.07.005, 2013.

%Vico, G., and Porporato, A.: Modelling C3 and C4 photosynthesis under water-stressed conditions, Plant Soil, 313,
%187-203, https://doi.org/10.1007/s11104-008-9691-4, 2008.


% RATIONALE:
% The following routine provides the stomatal conductance based on the co-limitation optimum model
% of Vico 2013 Ag Forest Met, but here revised to include the day respiration and leaf boundary layer
% conductances.
% The effects of declining soil and leaf water potential is accounted for via their effects on
% - the marginal water use efficiency lambda (influenced by psi_soil and not psi_leaf, because lambda does not
% respond istantaneously to changes in water availability, rather it responds more slowly), following Zhou Duursma et
% al 2013 Ag Forest Met and Manzoni et al 2011 Functional Ecol
% - the effects of leaf water potential on Vcmax and Jmax (following Vico and Porporato 2008 Plant and Soil).
% - the effects of leaf temperature on air humidity inside the leaf, so that, instead of VPD we use
% esat(Ts)-ea, and assume that surface temperature equals leaf temperature. In line with the use of Penman Monteith
% for the calculation of leaf temperature, using the approximation proposed by Penman, such difference is equated
% to esat(Ta)-eair+s(Tl-Ta)=VPD+s(Tl-Ta), where s is the slope of the curve
% relating saturation vapor pressure to temperature, assumed constant over the range Tl-Ta
%NB: there is an alternative version of this code, which includes just VPD, i.e., it assumes
%esat(Ts)=esat(Ta), or, in other words, that Tl=Ta. It is called  Opti_C3_RdgBL_psilVPD
%
% This stomatal model is coupled with a SPAC model, to determine (iteratively) the leaf water potential; this
% SPAC model includes the effects of soil and leaf water potentials on soil to root to leaf conductances,
% For this model, it is assumed that all the leaves are seing the same soil water content. This is the only
% way in which LAI is used in this code (and, as such, it does not matter how the leaves are organized in
% height, just their total amount)
%
% Notes on the optimization:
% 1) the leaf water potential and temperature are NOT included in the optimization directly, but they do affect Vcmax, Jmax and Lambda, i.e.,
% they are included in the whole model; this is the same approah Katul and Launiainen followed when
% implementing the optimization in APES. This is equivalent to assuming that the marginal effect of gs on Tl
% and psil is small with respect to those of gs on assimilation (but these assumptions apply only to the
% determination of stomatal conductance - see Manzoni et al 2011 Ecol Model for the inclusion of Tl and psil in the
% optimization, but without changes in Lambda)
% 2) cuticular conductance is neglected, on the ground that it becomes important only under extreme water stress
% (unrealistic for crops, which are anyhow managed so as not to reach extreme conditions)
% 3) no mesophyll conductance and its response to water availability is included - this would be impossible to
% separate from biochemical responses (see Zhou 2014 for a discussion on this aspect)
% 4) other responses (e.g., via changes in ABA - see Dewar) are also neglected

% Plant (and potentially species-specific) parameters:
% Vcmax25, Jmax25 - in input
% response of Vcmax and Jmax to water stress and Rd - default values in metabolic_waterstress
% Lambda under well watered conditions (in input or default in Lambda_waterstress), and its response to
% decreasing predawn leaf water potential (default response function and values in Lambda_waterstress, based
% on Manzoni et al 2011 and some recalculations; see Lambda_waterstress for details)
% gbl - in input
% Zr, RAI, LAI - in input
% dr - root diameter (default in soil2leafconductance)
% soil to root conductance, and its response to lowering water availability (two formulations possible - see soil2leaf conductance)
% root to leaf conductance, and its response to lowering water availability (two formulations possible - see soil2leaf conductance)

%Environmental parameters
% leaf temperature - in input
% PAR - in input
% VPD - in input

%Soil parameters
% pissat - in input
% Ksat - in input
% b - in input


% INPUTS:
%   Vcmax25: maximum carboxilation rate at 25 deg C and well watered cond (mumol/(m2leaf s))
%   Jmax25: maximum electron transport rate at 25 deg C and well watered cond (mumol/(m2leaf s))
%   Lambda_wwstar - lambda (mumol/mol): Lambda value under well watered conditions and at ca=400 ppm;
%       if negative, the value of Manzoni 2011 Func Ecol for forbes and herbs is used
%   T - leaf temperature (degC)
%   psilguess - vector containing two guesses for the leaf water potential, bracketing the solutuion - MPa, with sign
%   psipd - predawn leaf water potential, i.e., soil water potential at the beginning of the day (MPa) - with sign
%   psis - soil water potential (MPa) - with sign
%   Qp - PAR (mumol m-2 s-1)
%   ca - ambient CO2 concentration (ppm)
%   VPD - vapor pressure deficit (mol water/mol air)
%   Ta - air temperature (in deg C); if necessary, linkd VPD and Ta, but outside the function
%   gbl - leaf boundary layer conductance to water vapor (mol air/(m^2 leaf s))
%   ga - atmospheric conductance to water vapor or any other scalar (mol air/(m^2 leaf s))
%   psissat - soil water potential at saturation (MPa) - with sign
%   Ksat - soil hydraulic conductivity at saturation (m/s)
%   b - exponent of the soil water retention curve (-)
%   Zr - active rooting deoth (m)
%   RAI - root area index (-)
%   LAI - leaf area index (m2leaf/m2ground)
%   diagnostic: if 1, plots are made depicting the convergence of the numerical solution

% OUTPUTS:
%   An - net CO2 flux (mumol/(m2leaf s))
%   Eleaf - leaf transpiration rate (molH20/ (m2leaf s))
%   Eplant - transpiration rate per unit ground area (m3H20/(m2ground s))
%   gscopt - stomatal and cuticula conductances for vapor (mol air/(m2leaf s))
%   ci - leaf internal CO2 (ppm)
%   gsrp -

% FUNCTIONS CALLED:
% metabolic_waterstress
% Lambda_waterstress
% gsAciopt_RdgBL
% soil2leaf_conductance
% PARAMETERS
% CONSTANTS


PARAMETERS;
CONSTANTS;

%slope of the saturation vapour pressure-temperature relationship, calculated at the intermediate temp between
%Ta and Tl
Tinterm=(Tl+Ta)/2;
ss=4098.*(0.6108.*exp(17.27.*Tinterm./(Tinterm+237.3)))./(Tinterm+237.3).^2; %kPa/degC (NB This formula is slightly different from that reported in Jones, Appendix 4, but the numerical results are by and large the same, when considering this is in kPa/deg C
%humidity gradient, driving the transpiration
deltae=VPD+ss/Pair*(Tl-Ta); %P is the air pressure in kPa, so that deltae is in mol/mol, NOT in Pa (or kPa)



%Lambda responds to predawn water potential, i.e., psi soil, because it responds to slow(er) changes of leaf
%water potential than those occurring during the day; as such, Lambda_ws can be calculated only once per
%day
Lambda_ws=Lambda_waterstress(Lambda_wwstar, ca, psipd);

%In principle, if Ta>Tl and VPD low, deltae can be <0; in reality, this does not make sense, as high Ta means
%also high VPD; but a check on the sign of deltae is included anyhow and extremely low deltae are removed, or
%else complex results are obtained
if deltae>0.002 %minimum driving force of transpiration (below that, numerical issues occur)
    
    %leaf boundary layer and atmospheric conductances:
    gaCO2=ga; %atmospheric boundary layer conductance to CO2 is the same as that for vapor - 
              %this is based on turbulent transport, so it is the same for all scalars (including heat)
    gblCO2=gbl/bb; %leaf boundary layer conductance for CO2 needs a coeeficient, the ratio of diffusivities of vapor to CO2
    gblaCO2=gaCO2*gblCO2/(gaCO2+gblCO2);
    gbla=ga*gbl/(ga+gbl);
    unitconv=18*10^(-3)*10^(-3);
    %Notes on unit conversion:
    %gsla*deltae is in molH2O/(m^2leaf s)
    %gsrp is originally in  m^3H20/(m^2ground s MPa)
    %the unit conversion allows going from molH20 to kg (the 18*10^(-3)) and then from kgH20 to m^3liquidH20
    %(the second 10^(-3)); LAI added elsewhere will allow going from m^2leaf to m^2 ground
    
    %PLOTTING THE SUPPLY AND DEMAND CURVES - if run in diagnostic mode
    if diagnostic==1
        psils=[max(psilguess):-0.2:min(psilguess)];
        Eplants1=psils*0;Eplants2=psils*0;
        for i=1:length(psils)
            psil=psils(i);
            [Vcmax_ws,J_ws,Rd_ws,gamma_star,k1_ws,k2_ws]=metabolic_waterstress(Vcmax25, Jmax25,psil,Qp,Tl);
            %units: k1_c umolm-2s-1, k2_c: umol mol-1, gammast_c = umol mol-1
            
            %compute the optimal stomatal conductance
            [gsCO2opt,Aopt,ciopt]=gsAciopt_RdgBL(k1_ws,k2_ws,gamma_star,Rd_ws,Lambda_ws,ca,deltae,aa,gblaCO2);
            
            gcut=max(-gcutmax/psilstar*psil+gcutmax,0);
            gscopt=aa*gsCO2opt+gcut;
            gsla=gbla*gscopt/(gscopt+gbla); %stomatal, leaf boundary layer and bulk atmosphere conductance for vapor %mol air/(m^2 leaf s) - for layered model, gbla
            
            [gsrp,gsr,gp,gpalt]=soil2leaf_conductance(psil,psis,psissat,Ksat,bsoil,Zr,RAI);  %m/(s MPa) i.e. m^3H20/(m^2ground s MPa)
            Eleaf=gsla*deltae; %molH20/(m2leaf s);
            %in some formulas there is also a 0.622=18/29, i.e. ratio of the molecular weight of water and the efefctive one for dry air: we do not need this here because
            Eplants1(i)=Eleaf*unitconv; %m3H2O/(m2leaf s)
            Eplants2(i)=gsrp/LAI*(psis-psil); %m3H2O/(m2leaf s) - this is equivalent to assuming that all the leaves see the same soil to leaf conductacne
            
            %for use in the plot below - just to have a sense of the curves; assumes that psil=psis
            [gsrps(i),gsrs(i),gps(i),gpsalt(i)]=soil2leaf_conductance(psils(i),psils(i),psissat,Ksat,bsoil,Zr,RAI);
            Lambda_wss(i)=Lambda_waterstress(Lambda_wwstar, ca, psils(i));
            Vcmax_wss(i)=Vcmax_ws;
            J_wss(i)=J_ws;
        end
        figure(101)
        plot(psils,real(Eplants1),'-g');hold on
        plot(psils,imag(Eplants1),':g');
        plot(psils,real(Eplants2),'-r');
        plot(psils,imag(Eplants2),':r');
        title('green: demand, red: supply; solid: real part, dotted: imaginary')
        ylabel('E (m^3_{H2O}/(m^2_{leaf} s)')
        xlabel('\psi_L (MPa)')
        
        figure(102)
        subplot(2,3,1)
        plot(psils,Vcmax_wss,'-')
        xlabel('\psi_L (MPa)');ylabel('V_{c,max} (\mu mol/(m^2_{leaf} s))')
        subplot(2,3,2)
        plot(psils,J_wss,'-')
        xlabel('\psi_L (MPa)');ylabel('J (\mu mol/(m^2_{leaf} s))')
        subplot(2,3,3)
        plot(psils,gps,'-');hold on
        plot(psils,gpsalt,':');
        xlabel('\psi_L (MPa)');ylabel('g_P (m/(s MPa))')
        subplot(2,3,4)
        plot(psils,gsrs,'-');
        xlabel('\psi_S (MPa)');ylabel('g_{SR} (m/(s MPa))')
        subplot(2,3,5)
        Lambda_ws=Lambda_waterstress(Lambda_wwstar, ca, psipd);
        semilogy(psipd,Lambda_ws,'*');hold on
        semilogy(psils,Lambda_wss,'-')
        xlabel('\psi_S (MPa)');ylabel('\lambda (mol/mol)')
        hold off
        
    end
    
    clear Vcmax_ws J_ws Rd_ws gamma_star k1_ws k2_ws gsCO2opt Aopt ciopt
    
    %NUMERICAL SOLUTION
    %based on the bisection method, to ensure convergence
    %the routine is exited when the change in the resulting psil is smaller than tolerPSIL, or after a number of
    %iterations equal to maxiter
    iter=1;  
    
    % initial guesses for the bisection method
    psil1=min(psilguess);psil2=min(max(psilguess),psis);
    psil3=(psil1+psil2)/2;
    f3=1;   
    
    while (abs(real(psil2-psil1))>tolerpsil && iter<maxiter) || iter<2
        psil=psil1;
        %compute metabolic paramters and lambda, including the effects of water stress and leaf temperature
        [Vcmax_ws,J_ws,Rd_ws,gamma_star,k1_ws,k2_ws]=metabolic_waterstress(Vcmax25, Jmax25,psil,Qp,Tl); %units: k1_ws: mumol m-2s-1, k2_ws: mumol mol-1, gammast_star: mumol mol-1
        
        %compute the optimal stomatal conductance
        [gsCO2opt,Aopt,ciopt]=gsAciopt_RdgBL(k1_ws,k2_ws,gamma_star,Rd_ws,Lambda_ws,ca,deltae,aa,gblaCO2);
        
        gcut=max(-gcutmax/psilstar*psil+gcutmax,0);
        gscopt=aa*gsCO2opt+gcut;
        gsla=gbla*gscopt/(gscopt+gbla); %stomatal to atmosphere conductance for vapor %mol air/(m^2 leaf s)
        [gsrp,gsr,gp,gpalt]=soil2leaf_conductance(psil,psis,psissat,Ksat,bsoil,Zr,RAI);  %m/(s MPa) i.e. m^3H20/(m^2ground s MPa)
        Eplant2=gsrp/LAI*(psis-psil); %it is assumed that all the leaves see the same soil to leaf conductance
        Eleaf=gsla.*deltae; Eplant1=Eleaf*unitconv;
        f1=Eplant1-Eplant2;Eplant11=Eplant1;Eplant21=Eplant2;
        
        psil=psil2;
        %compute metabolic paramters and lambda, including the effects of water stress and leaf temperature
        [Vcmax_ws,J_ws,Rd_ws,gamma_star,k1_ws,k2_ws]=metabolic_waterstress(Vcmax25, Jmax25,psil,Qp,Tl); %units: k1_ws: mumol m-2s-1, k2_ws: mumol mol-1, gammast_star: mumol mol-1
        %compute the optimal stomatal conductance
        [gsCO2opt,Aopt,ciopt]=gsAciopt_RdgBL(k1_ws,k2_ws,gamma_star,Rd_ws,Lambda_ws,ca,deltae,aa,gblaCO2);
        gcut=max(-gcutmax/psilstar*psil+gcutmax,0);
        gscopt=aa*gsCO2opt+gcut;
        gsla=gbla*gscopt/(gscopt+gbla); %stomatal to atmosphere conductance for vapor %mol air/(m^2 leaf s)
        [gsrp,gsr,gp,gpalt]=soil2leaf_conductance(psil,psis,psissat,Ksat,bsoil,Zr,RAI);  %m/(s MPa) i.e. m^3H20/(m^2ground s MPa)
        Eplant2=gsrp/LAI*(psis-psil); %it is assumed that all the leaves see the same soil to leaf conductance
        Eleaf=gsla.*deltae; Eplant1=Eleaf*unitconv;
        f2=Eplant1-Eplant2;Eplant12=Eplant1;Eplant22=Eplant2;
        
        if iter==1
            if f1*f2>0
                psil1
                f1
                Eplant11
                Eplant21
                psil2
                f2
                Eplant12
                Eplant22
                psis
                pause
                %if the first set of guesses does not work, then it tries with another one, before giving up;
                %in all cases, psil<psis, for physical reasons; yet, posing psil2=psis may not work in some cases,
                %so both psil2=psis and psis=0 are explored as upper guess
                if max(psilguess)<0
                    psil2=0;
                else
                    psil2=psis-0.1;
                end
                psil=psil1;
                %compute metabolic paramters and lambda, including the effects of water stress and leaf temperature
                [Vcmax_ws,J_ws,Rd_ws,gamma_star,k1_ws,k2_ws]=metabolic_waterstress(Vcmax25, Jmax25,psil,Qp,Tl); %units: k1_ws: mumol m-2s-1, k2_ws: mumol mol-1, gammast_star: mumol mol-1
                %compute the optimal stomatal conductance
                [gsCO2opt,Aopt,ciopt]=gsAciopt_RdgBL(k1_ws,k2_ws,gamma_star,Rd_ws,Lambda_ws,ca,deltae,aa,gblaCO2);
                gcut=max(-gcutmax/psilstar*psil+gcutmax,0);
                gscopt=aa*gsCO2opt+gcut;
                gsla=gbla*gscopt/(gscopt+gbla); %stomatal to atmosphere conductance for vapor %mol air/(m^2 leaf s)
                [gsrp,gsr,gp,gpalt]=soil2leaf_conductance(psil,psis,psissat,Ksat,bsoil,Zr,RAI);  %m/(s MPa) i.e. m^3H20/(m^2ground s MPa)
                Eplant2=gsrp/LAI*(psis-psil); %it is assumed that all the leaves see the same soil to leaf conductance
                Eleaf=gsla.*deltae; Eplant1=Eleaf*unitconv;
                f1=Eplant1-Eplant2;Eplant11=Eplant1;Eplant21=Eplant2;
                
                psil=psil2;
                %compute metabolic parameters and lambda, including the effects of water stress and leaf temperature
                [Vcmax_ws,J_ws,Rd_ws,gamma_star,k1_ws,k2_ws]=metabolic_waterstress(Vcmax25, Jmax25,psil,Qp,Tl); %units: k1_ws: mumol m-2s-1, k2_ws: mumol mol-1, gammast_star: mumol mol-1
                %compute the optimal stomatal conductance
                [gsCO2opt,Aopt,ciopt]=gsAciopt_RdgBL(k1_ws,k2_ws,gamma_star,Rd_ws,Lambda_ws,ca,deltae,aa,gblaCO2);
                gcut=max(-gcutmax/psilstar*psil+gcutmax,0);
                gscopt=aa*gsCO2opt+gcut;
                gsla=gbla*gscopt/(gscopt+gbla); %stomatal to atmosphere conductance for vapor %mol air/(m^2 leaf s)
                [gsrp,gsr,gp,gpalt]=soil2leaf_conductance(psil,psis,psissat,Ksat,bsoil,Zr,RAI);  %m/(s MPa) i.e. m^3H20/(m^2ground s MPa)
                Eplant2=gsrp/LAI*(psis-psil); %it is assumed that all the leaves see the same soil to leaf conductance
                Eleaf=gsla.*deltae; Eplant1=Eleaf*unitconv;
                f2=Eplant1-Eplant2;Eplant12=Eplant1;Eplant22=Eplant2;
                if f1*f2>0 %if even the second round of guesses does not work, than it throws in the towel
                    psil1
                    f1
                    Eplant11
                    Eplant21
                    psil2
                    f2
                    Eplant12
                    Eplant22
                    psis
                    display('initial guesses not appropriate')
                    return
                end                
            end
        end
        
        psil=(psil1+psil2)/2;psil3=psil;
        %compute metabolic paramters and lambda, including the effects of water stress and leaf temperature
        [Vcmax_ws,J_ws,Rd_ws,gamma_star,k1_ws,k2_ws]=metabolic_waterstress(Vcmax25, Jmax25,psil,Qp,Tl);   %units: k1_ws: mumol m-2s-1, k2_ws: mumol mol-1, gammast_star: mumol mol-1
        %compute the optimal stomatal conductance
        [gsCO2opt,Aopt,ciopt]=gsAciopt_RdgBL(k1_ws,k2_ws,gamma_star,Rd_ws,Lambda_ws,ca,deltae,aa,gblaCO2);
        gcut=max(-gcutmax/psilstar*psil+gcutmax,0);
        gscopt=aa*gsCO2opt+gcut;
        gsla=gbla*gscopt/(gscopt+gbla); %stomatal to atmosphere conductance for vapor %mol air/(m^2 leaf s)
        [gsrp,gsr,gp,gpalt]=soil2leaf_conductance(psil,psis,psissat,Ksat,bsoil,Zr,RAI);  %m/(s MPa), i.e. m^3H20/(m^2ground s MPa)
        Eplant2=gsrp/LAI*(psis-psil); %it is assumed that all the leaves see the same soil to leaf conductance
        Eleaf=gsla.*deltae; Eplant1=Eleaf*unitconv;
        
        f3=Eplant1-Eplant2;Eplant13=Eplant1;Eplant23=Eplant2;
        
        if diagnostic==1
            figure(101)
            plot(psil1,real(Eplant21),'o');
            plot(psil1,real(Eplant11),'o');
            plot(psil2,real(Eplant22),'v')
            plot(psil2,real(Eplant12),'v')
            plot(psil3,real(Eplant23),'*')
            plot(psil1,real(Eplant13),'*')
        end
        
        if f1*f3<0
            psil1=psil1;
            psil2=psil3;
        else
            psil1=psil3;
            psil2=psil2;
        end
        iter=iter+1;
    end
    
    if diagnostic==1
        if iter<maxiter
            display('iterations needed for convergence to be reached')
            iter
        else
            display('convergence was not reached')
        end
    end
    
    
    psil=psil3;
    if isfinite(ciopt)==1
        Anleaf=Aopt; %molCO2/(m2leaf s)
    else
        Anleaf=-Rd_ws;
    end
    Eleaf=gsla.*deltae; %molH20/(m2leaf s)
    
    Eplant2=gsrp/LAI*(psis-psil); %it is assumed that all the leaves see the same soil to leaf conductance  %m3H2O/(m2leaf s)
    Eleaf=gsla.*deltae; Eplant1=Eleaf*unitconv; %m3H2O/(m2leaf s)
    
    if abs(real(Eplant1-Eplant2))>10^(-8)
        display('something is amiss')
        psis
        iter
        gsrp
        real(Eplant2)
        Eplant=Eplant1
        %pause
    else
        Eplant=real(Eplant1);
    end
    
    if diagnostic==1
        figure(101)
        plot(psil,Eplant,'hr')
        hold off
    end
    
else
    if diagnostic==1
        deltae
        display('all results set to 0')
    end
    Anleaf=0; Eleaf=0; Eplant=0; gsCO2opt=0; ciopt=0; psil=psis;Vcmax_ws=Vcmax25;
end

%Uncomment the following line for debugging purposes
%gsCO2opt

return

