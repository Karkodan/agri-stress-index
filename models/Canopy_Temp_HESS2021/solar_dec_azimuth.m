function [ZEN, AZIM,dd]=solar_dec_azimuth (DOY, LAT1, tm)

%Author: Giulia Vico
%Date: September 2018
%This is part of the codes used to generate the results in the following pubblication.
%Please refer to that for details on the model, its rationale and parameterization.
%Luan, X. and Vico, G. (2021), Canopy temperature and heat stress are increased by compound high air temperature
%and water stress, and reduced by irrigation ? A modeling analysis, Hydrol. Earth Syst. Sci.,
%https://doi.org/10.5194/hess-2020-549

%Key reference: 
%Dingman, S.: Physical hydrology, Macmillan, New York, USA, 575 pp., 1994.

%INPUTS
%DOY day of the year of interest
%LAT latitude (in deg)
%tm solar time of interest (hrs)

%OUTPUTS
%ZEN zenith angle
%AZIM azimuth
%dd solar declination angle

%FUNCTIONS CALLED
%none

%	----> compute Solar declination angle  
CF=pi/180;
LAT=LAT1*CF;
xx=278.97+0.9856*DOY+1.9165*sin((356.6+0.9856*DOY)*CF);
dd=asin(0.39785*sin(xx*CF));
%	----> compute Zenith angle
f=(279.575+0.9856*DOY)*CF;
ET=(-104.7*sin(f)+596.2*sin(2*f)+4.3*sin(3*f)-12.7*sin(4*f)-429.3*cos(f)-2*cos(2*f)+19.3*cos(3*f))/3600;
    
TF=ET;
aa=sin(LAT)*sin(dd)+(cos(LAT))*cos(dd)*cos(15*(tm-12-TF)*CF);
ZEN=acos(aa);
%	----> compute azimuth angle	
AZIM=acos(-(sin(dd)-(sin(LAT))*cos(ZEN))/((cos(LAT))*sin(ZEN)));

return
	
