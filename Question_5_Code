%% Load Data
load('Data.mat')

% Assuming the length of your data is even
half_length = length(D) / 2;

% Calibration data
D_calib = D(1:half_length);
P_calib = P(1:half_length);
Q_calib = Q(1:half_length);
T_calib = T(1:half_length);
Tmax_calib = Tmax(1:half_length);
Tmin_calib = Tmin(1:half_length);

% Validation data
D_valid = D(half_length+1:end);
P_valid = P(half_length+1:end);
Q_valid = Q(half_length+1:end);
T_valid = T(half_length+1:end);
Tmax_valid = Tmax(half_length+1:end);
Tmin_valid = Tmin(half_length+1:end);

%% Catchment properties
Lat = 40; % Latitude [deg]

% Calibration discharge for comparison
Qobs_calib = Q_calib;

%% Model parameters
T_T = 0; % Temperature threshold [deg C] distinguishing between rain and snow
Tsm = 0; % Temperature threshold for melting [deg C]
k = 2; % Daily degree-day snowmelt parameter [mm/day/degC]

%% Calibration: Partition snow and rainfall
SN_calib = zeros(length(D_calib),1); 
R_calib = zeros(length(D_calib),1);

SN_calib(Tmax_calib<=T_T) = P_calib(Tmax_calib<=T_T);
R_calib(Tmin_calib>T_T) = P_calib(Tmin_calib>T_T);
R_calib(Tmin_calib<=T_T & Tmax_calib>T_T) = ...
    (Tmax_calib(Tmin_calib<=T_T & Tmax_calib>T_T).*P_calib(Tmin_calib<=T_T & Tmax_calib>T_T)) ./ ...
    (Tmax_calib(Tmin_calib<=T_T & Tmax_calib>T_T)-Tmin_calib(Tmin_calib<=T_T & Tmax_calib>T_T));
SN_calib(Tmin_calib<=T_T & Tmax_calib>T_T) = P_calib(Tmin_calib<=T_T & Tmax_calib>T_T) - R_calib(Tmin_calib<=T_T & Tmax_calib>T_T);
SN_calib(SN_calib<0) = 0;

%% Calibration: Snow degree day model
SC_calib = zeros(length(D_calib),1); 
M_calib = zeros(length(D_calib),1);

for t = 2:length(D_calib)
    SC_calib(t) = SC_calib(t-1) + SN_calib(t);
    if T_calib(t) > Tsm
        M_calib(t) = min(k*(T_calib(t)-Tsm),SC_calib(t));
    else
        M_calib(t) = 0;
    end
    SC_calib(t) = SC_calib(t) - M_calib(t);
end

%% Calibration: Compute PET with routine by Hamon
Jday_calib = day(datetime(datestr(D_calib)),'dayofyear');
PET_calib = zeros(length(D_calib),1);

for t = 1:length(D_calib)
    delta = 0.4093*sin((2*pi/365)*Jday_calib(t)-1.405);
    omega_s = acos(-tan(2*pi*Lat/360).*tan(delta));
    Nt = 24*omega_s/pi;
    es = 0.6108*exp(17.27*T_calib(t)./(T_calib(t)+237.3));
    PET_calib(t) =(2.1*(Nt.^2).*es)./(T_calib(t)+273.3);
end

%% Calibration: Soil water content
Smax_vals = 20:5:100; % Values for soil storage capacity
fr_vals = 0.01:0.03:0.2; % Values for runoff fraction 
rs_vals = 5:5:100; % Values for soil reservoir constant
rg_vals = 5:5:100; % Values for groundwater reservoir constant

Ind = 1;
for Smax = Smax_vals
    for fr = fr_vals
        for rs = rs_vals
            for rg = rg_vals
                AET_calib = zeros(length(D_calib),1); 
                S_calib = zeros(length(D_calib),1); 
                Qs_calib = zeros(length(D_calib),1);
                I_calib = zeros(length(D_calib),1);

                for t = 2:length(D_calib)
                    if S_calib(t-1)+R_calib(t)+M_calib(t) > Smax
                        Qs_calib(t) = (R_calib(t)+M_calib(t))*fr;
                    end
                    S_calib(t) = S_calib(t-1) + M_calib(t) + R_calib(t) - Qs_calib(t);

                    AET_calib(t) = (S_calib(t)/Smax)*PET_calib(t);
                    if S_calib(t) > Smax
                        AET_calib(t) = PET_calib(t);
                    end

                    S_calib(t) = S_calib(t) - AET_calib(t);

                    I_calib(t) = (1/rs)*S_calib(t);
                    if S_calib(t) > Smax
                        I_calib(t) = S_calib(t) - Smax;
                    end
                    S_calib(t) = S_calib(t) - I_calib(t);
                end

                %% Calibration: Groundwater and streamflow
                GW_calib = zeros(length(D_calib),1); 
                Qgw_calib = zeros(length(D_calib),1); 
                Qsim_calib = zeros(length(D_calib),1);

                for t = 2:length(D_calib)
                    GW_calib(t) = GW_calib(t-1) + I_calib(t);
                    Qgw_calib(t) = (1/rg)*GW_calib(t);
                    GW_calib(t) = GW_calib(t) - Qgw_calib(t);
                    Qsim_calib(t) = Qgw_calib(t) + Qs_calib(t);
                end

                %% Calibration: Evaluation Metrics (KGE)
                KGE(Ind,1) = 1 - sqrt((corr(Qsim_calib,Qobs_calib)-1)^2 + ...
                    (std(Qsim_calib)/std(Qobs_calib)-1)^2 + ...
                    (mean(Qsim_calib)/mean(Qobs_calib)-1)^2);

                Par(Ind,:) = [Smax fr rs rg]; 
                Ind = Ind + 1;
            end
        end
    end
end

%% Calibration: Finding the Best Parameter Sets
[bestKGE, bestKGEind] = max(KGE);
bestParKGE = Par(bestKGEind, :);

%% Validation: Partition snow and rainfall
SN_valid = zeros(length(D_valid),1); 
R_valid = zeros(length(D_valid),1);

SN_valid(Tmax_valid<=T_T) = P_valid(Tmax_valid<=T_T);
R_valid(Tmin_valid>T_T) = P_valid(Tmin_valid>T_T);
R_valid(Tmin_valid<=T_T & Tmax_valid>T_T) = ...
    (Tmax_valid(Tmin_valid<=T_T & Tmax_valid>T_T).*P_valid(Tmin_valid<=T_T & Tmax_valid>T_T)) ./ ...
    (Tmax_valid(Tmin_valid<=T_T & Tmax_valid>T_T)-Tmin_valid(Tmin_valid<=T_T & Tmax_valid>T_T));
SN_valid(Tmin_valid<=T_T & Tmax_valid>T_T) = P_valid(Tmin_valid<=T_T & Tmax_valid>T_T) - R_valid(Tmin_valid<=T_T & Tmax_valid>T_T);
SN_valid(SN_valid<0) = 0;

%% Validation: Snow degree day model
SC_valid = zeros(length(D_valid),1); 
M_valid = zeros(length(D_valid),1);

for t = 2:length(D_valid)
    SC_valid(t) = SC_valid(t-1) + SN_valid(t);
    if T_valid(t) > Tsm
        M_valid(t) = min(k*(T_valid(t)-Tsm),SC_valid(t));
    else
        M_valid(t) = 0;
    end
    SC_valid(t) = SC_valid(t) - M_valid(t);
end

%% Validation: Compute PET with routine by Hamon
Jday_valid = day(datetime(datestr(D_valid)),'dayofyear');
PET_valid = zeros(length(D_valid),1);

for t = 1:length(D_valid)
    delta = 0.4093*sin((2*pi/365)*Jday_valid(t)-1.405);
    omega_s = acos(-tan(2*pi*Lat/360).*tan(delta));
    Nt = 24*omega_s/pi;
    es = 0.6108*exp(17.27*T_valid(t)./(T_valid(t)+237.3));
    PET_valid(t) = (2.1*(Nt.^2).*es)./(T_valid(t)+273.3);
end

%% Validation: Soil water content
AET_valid = zeros(length(D_valid),1); 
S_valid = zeros(length(D_valid),1); 
Qs_valid = zeros(length(D_valid),1);
I_valid = zeros(length(D_valid),1);

for t = 2:length(D_valid)
    if S_valid(t-1)+R_valid(t)+M_valid(t) > Smax
        Qs_valid(t) = (R_valid(t)+M_valid(t))*fr;
    end
    S_valid(t) = S_valid(t-1) + M_valid(t) + R_valid(t) - Qs_valid(t);

    AET_valid(t) = (S_valid(t)/Smax)*PET_valid(t);
    if S_valid(t) > Smax
        AET_valid(t) = PET_valid(t);
    end

    S_valid(t) = S_valid(t) - AET_valid(t);

    I_valid(t) = (1/rs)*S_valid(t);
    if S_valid(t) > Smax
        I_valid(t) = S_valid(t) - Smax;
    end
    S_valid(t) = S_valid(t) - I_valid(t);
end

%% Validation: Groundwater and streamflow
GW_valid = zeros(length(D_valid),1); 
Qgw_valid = zeros(length(D_valid),1); 
Qsim_valid = zeros(length(D_valid),1);

for t = 2:length(D_valid)
    GW_valid(t) = GW_valid(t-1) + I_valid(t);
    Qgw_valid(t) = (1/rg)*GW_valid(t);
    GW_valid(t) = GW_valid(t) - Qgw_valid(t);
    Qsim_valid(t) = Qgw_valid(t) + Qs_valid(t);
end

%% Validation: Evaluation Metrics (NSE)
Qobs_valid = Q_valid;
NSE_valid = 1 - ( sum( (Qobs_valid - Qsim_valid).^2 ) / sum( (Qobs_valid - mean(Qobs_valid)).^2 ) );

%% Validation: Finding the Best Parameter Sets
bestNSE = NSE_valid;
bestParNSE = Par(bestKGEind, :);

% Display results
disp("Best parameter set using KGE: " + num2str(bestParKGE))
disp("Best parameter set using NSE: " + num2str(bestParNSE))
