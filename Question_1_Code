% WatBal model - Watershed modeling course, UNIL 2024
%% Input
% D - Dates, Matlab format
% P - Daily precipitation [mm]
% Q - Daily discharge [mm/day]
% T, Tmax and Tmin - Mean, maximum, and minimum daily temperature [deg Celsius]
% Jday - Day of the year (Julian date)
load('Data.mat')
%% Catchment properties
Lat = 40; % Latitude [deg]
Qobs = Q; % Observed discharge [mm/day]
%% Model parameters
T_T = 0; % Temperature threshols [deg C] distingishing between rain and snow
Tsm = 0; % Temperature threshols for melting [deg C]
k = 2; % Daily degree-day snowmelt parameter [mm/day/degC]
Smax = 50; % Soil storage capacity, maximum soil water content [mm]
fr = 0.15; % Fraction of watershed area that is generating surface runoff [0-1];
rs = 35; % Soil water linear reservoir constant, mean residence time [days]
rg = 35; % Ground water linear reservoir constant, mean residence time [days]
%% Partition snow and rainfall:
SN=zeros(length(D),1); R=zeros(length(D),1); % Creating snow (SN) and rainfall (R) time series with equal length of D
SN(Tmax<=T_T)=P(Tmax<=T_T);
R(Tmin>T_T)=P(Tmin>T_T);
R(Tmin<=T_T&Tmax>T_T)=(Tmax(Tmin<=T_T&Tmax>T_T).*P(Tmin<=T_T&Tmax>T_T))./(Tmax(Tmin<=T_T&Tmax>T_T)-Tmin(Tmin<=T_T&Tmax>T_T));
SN(Tmin<=T_T&Tmax>T_T)=P(Tmin<=T_T&Tmax>T_T)-R(Tmin<=T_T&Tmax>T_T); SN(SN<0)=0;
%% Snow degree day model
SC=zeros(length(D),1); M=zeros(length(D),1);
for t = 2:length(D)
    % Compute snow cover (SC)
    SC(t) = SC(t-1) + SN(t);
    % Compute snow melt
    if T(t) > Tsm
        M(t) = min(k*(T(t)-Tsm),SC(t));
    else
        M(t) = 0;
    end
    % Compute snow cover after melt
    SC(t) = SC(t) - M(t);
end
%% Compute PET with routine by Hamon
% Create a variable named 'Jday' which contains the dates in Julian days
Jday=day(datetime(datestr(D)),'dayofyear');
PET=zeros(length(D),1);
for t = 1:length(D)
    delta = 0.4093*sin((2*pi/365)*Jday(t)-1.405);
    omega_s = acos(-tan(2*pi*Lat/360).*tan(delta));
    Nt = 24*omega_s/pi;
    es = 0.6108*exp(17.27*T(t)./(T(t)+237.3));
    PET(t) =(2.1*(Nt.^2).*es)./(T(t)+273.3);
end
%% Soil water content
AET=zeros(length(D),1); S=zeros(length(D),1); Qs=zeros(length(D),1);
I=zeros(length(D),1);
for t = 2:length(D)
    % Compute surface runoff
    if S(t-1)+R(t)+M(t) > Smax
        Qs(t) = (R(t)+M(t))*fr;
    end
    % Compute water balance
    S(t) = S(t-1) + M(t) + R(t) - Qs(t);
    % Compute AET
    AET(t)=(S(t)/Smax)*PET(t);
    if S(t) > Smax
        AET(t)=PET(t);
    end
    % Remove AET from the soil water content
    S(t) = S(t)-AET(t);
    % Compute infiltration to soil reservoir
    I(t)=(1/rs)*S(t);
    if S(t) > Smax
        I(t) = S(t)-Smax;
    end
    S(t)=S(t)-I(t);
end
%% Groundwater and streamflow
GW=zeros(length(D),1); Qgw=zeros(length(D),1); Qsim=zeros(length(D),1);
for t = 2:length(D)
    GW(t) = GW(t-1) + I(t);
    % Compute infiltration to groundwater reservoir
    Qgw(t)=(1/rg)*GW(t);
    GW(t)=GW(t)-Qgw(t);
    Qsim(t)=Qgw(t)+Qs(t);
end
%% KGE
KGE=1-sqrt((corr(Qsim,Qobs)-1)^2+(std(Qsim)/std(Qobs)-1)^2+(mean(Qsim)/mean(Qobs)-1)^2);
%% NSE
NSE = 1-sum((Qsim-Qobs).^2)/sum((Qobs-mean(Qobs)).^2);

