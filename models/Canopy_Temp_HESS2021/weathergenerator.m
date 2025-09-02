function [rain,T]=weathergenerator(dt,Tmaxmean,Tin,rainY,Fr,N,k1,k3)

%Author: Xiangyu Luan and Giulia Vico 
%Date: December 2020
%Part of the codes used to generate the results in the following pubblication. 
%Please refer to that for details on the model, its rationale and parameterization.
%Luan, X. and Vico, G. (2021), Canopy temperature and heat stress are increased by compound high air temperature 
%and water stress, and reduced by irrigation ? A modeling analysis, Hydrol. Earth Syst. Sci., 
%https://doi.org/10.5194/hess-2020-549 

%INPUTS:
% Tmean: the mean air temperature for simulation period
% Tin: intial air temperature (used in the first time step)
% rainY: total average annual precipitation (mm)
% Fr: average precipitation frequency (1/d)
% k1: mean reversion rate of the Ornstein Uhlenbeck process, to generate air temperature (i.e., 1/k1 is the relaxation time) (1/d)
% k3: diffusion parameter of the Ornstein Uhlenbeck process (deg C^2/d)

%OPUTPUT:
% rain: time series of the daily precipitations (mm)
% T: time series of the daily air temperature (deg C)


%FUNCTIONS CALLED:
%none


% ____________________________precipitation  generator__________________
nr=round(N*Fr);  %Total rainfall events
a=N/365;
alpha=rainY/nr*a;% alpha is the mean depth of rainfall events
lambda=1/Fr;
rain=zeros(1,N);
for i=1:N
    rain(i)=0;
end
pd_tau=makedist('Exponential','mu',lambda); % generates an exponential interarrival time, tau
pd_h=makedist('Exponential','mu',alpha);    % generates an exponetial rainfall depth
tau=random(pd_tau);      %tau=time of the new rainfall event

while tau<N
    rain_day=ceil(tau);
    rain(rain_day)=rain(rain_day)+random(pd_h);
    tau=tau+random(pd_tau);
end


% ____________________________Air temperature generator__________________
[T]=zeros(1,N);%preset T
T(1)=Tin;

for i=1:N-1
    dT=(-k1*(T(i)-Tmaxmean)+sqrt(k3)*randn(1,1))*dt;
    T(i+1)=T(i)+dT;
end




return