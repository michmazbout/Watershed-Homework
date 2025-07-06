% WatBal model - Watershed modeling course, UNIL 2024
% Revised version with parameter calibration

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

% Model parameters to be calibrated
k_range = 1:0.5:5; % Define a range with smaller increments for k
T_T_range = -5:1:5; % Define a range with smaller increments for T_T

% Initialize variables to store best parameter set and performance
best_KGE = -inf;
best_NSE = -inf;
best_params = [];
best_Qsim = [];

% Loop through parameter combinations for calibration
for k = k_range
    for T_T = T_T_range
        
        % Snow and rainfall partitioning
        SN=zeros(length(D),1); 
        R=zeros(length(D),1); 
        
        SN(Tmax<=T_T)=P(Tmax<=T_T);
        R(Tmin>T_T)=P(Tmin>T_T);
        R(Tmin<=T_T & Tmax>T_T)=(Tmax(Tmin<=T_T & Tmax>T_T).*P(Tmin<=T_T & Tmax>T_T))./(Tmax(Tmin<=T_T & Tmax>T_T)-Tmin(Tmin<=T_T & Tmax>T_T));
        SN(Tmin<=T_T & Tmax>T_T)=P(Tmin<=T_T & Tmax>T_T)-R(Tmin<=T_T & Tmax>T_T); 
        SN(SN<0)=0;
        
        % Snow degree day model
        SC=zeros(length(D),1); 
        M=zeros(length(D),1);
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
        
        % PET computation
        Jday=day(datetime(datestr(D)),'dayofyear');
        PET=zeros(length(D),1);
        for t = 1:length(D)
            delta = 0.4093*sin((2*pi/365)*Jday(t)-1.405);
            omega_s = acos(-tan(2*pi*Lat/360).*tan(delta));
            Nt = 24*omega_s/pi;
            es = 0.6108*exp(17.27*T(t)./(T(t)+237.3));
            PET(t) =(2.1*(Nt.^2).*es)./(T(t)+273.3);
        end
        
        % Soil water content
        Smax_vals = 20:5:100; % Values for soil storage capacity
        fr_vals = 0.01:0.03:0.2; % Values for runoff fraction 
        rs_vals = 5:5:100; % Values for soil reservoir constant
        rg_vals = 5:5:100; % Values for groundwater reservoir constant
        
        Ind = 1;
        for Smax = Smax_vals
            for fr = fr_vals
                for rs = rs_vals
                    for rg = rg_vals
                        AET=zeros(length(D),1); 
                        S=zeros(length(D),1); 
                        Qs=zeros(length(D),1);
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
                        
                        % Groundwater and streamflow
                        GW=zeros(length(D),1); 
                        Qgw=zeros(length(D),1); 
                        Qsim=zeros(length(D),1);
                        for t = 2:length(D)
                            GW(t) = GW(t-1) + I(t);
                            % Compute infiltration to groundwater reservoir
                            Qgw(t)=(1/rg)*GW(t);
                            GW(t)=GW(t)-Qgw(t);
                            Qsim(t)=Qgw(t)+Qs(t);
                        end
                        
                        % KGE and NSE calculation
                        KGE = 1-sqrt((corr(Qsim,Qobs)-1)^2+(std(Qsim)/std(Qobs)-1)^2+(mean(Qsim)/mean(Qobs)-1)^2);
                        NSE = 1 - ( sum( (Qobs - Qsim).^2 ) / sum( (Qobs - mean(Qobs)).^2 ) ); 
                        
                        % Update best parameter set if performance improves
                        if KGE > best_KGE && NSE > best_NSE
                            best_KGE = KGE;
                            best_NSE = NSE;
                            best_params = [k, T_T, Smax, fr, rs, rg];
                            best_Qsim = Qsim;
                        end
                        
                        Ind = Ind + 1;
                    end
                end
            end
        end
    end
end

%% Display and report the best parameter set achieved
disp('Best parameter set:')
disp(['k = ', num2str(best_params(1))])
disp(['T_T = ', num2str(best_params(2))])
disp(['Smax = ', num2str(best_params(3))])
disp(['fr = ', num2str(best_params(4))])
disp(['rs = ', num2str(best_params(5))])
disp(['rg = ', num2str(best_params(6))])

%% Plot observed vs. simulated streamflow
figure;
plot(D, Qobs, 'b', D, best_Qsim, 'r');
xlabel('Date');
ylabel('Streamflow (mm/day)');
title('Observed vs. Simulated Streamflow');
legend('Observed', 'Simulated');
