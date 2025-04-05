%clear
%clc
%cd 'C:\Users\rakif\Box\DBF 2023-2024\Sim Dev'
%% Plotting Settings
c = [0.647 0.078 0.090;
     0.424 0.451 0.451;
     0.784 0.784 0.784;
     0 0.451 0.376;
     0 0.373 0.522;
     1 0.8 0];
set(0,'DefaultAxesColorOrder',c)
set(0,'defaultLineLineWidth',1.5)
%%
% Read column B of Inputs excel sheet
variables = readmatrix('optimizerV2_inputs.xlsx','Range','B:B');
Cd_min = variables(1); % from excel
Cl_max = variables(2); % from excel
g = 9.81;
k = variables(3); % from excel
eff_c = 0.75;
eff_e = 0.9;
mu = 0.04;
rho = 1.225;
lap_dist = 762;
eta_coeff = [variables(4) variables(5) variables(6) variables(7)]; % from excel
n = 25; %% Do not go above 40
%Constraints
TOFL_max = 5.5;
WL_max = 98.1;
P_max = 2000;
S_max = 0.5;
passenger_max = 50;

% Ground Mission Timing
Time2change_config = 10; % s; Estimated time to change plane configuration from parking to flight and vice-versa
Time2load_m2_payload = 30; % s; Estimated time to load entire m2 payload
Time2unload_m2_payload = 30; % s; Estimated time to unload entire m2 paylaod, including crew (which is removed later)
Time2secure_hatches_doors = 5; % s; Estimated time to secure hatches and doors
Time2unsecure_hatches_doors = 5; % s; Estimated time to unsecure hatches and doors
% BEST_CONFIG DEPENDS ON THESE 2 VALUES BELOW:
Time2loadunload_passenger = 1;

Sweep_S = linspace(0.1, 1, n);
Sweep_W_p = linspace(0, 5, n);
Sweep_v_max = linspace(15,45,n);
Sweep_passengers = linspace(1,3*n-2,n);
Sweep_C = linspace(32,100,n)*3600;
% Coordinates go (S, W_p, v_max, passengers)
tic
[S, W_p, v_max, passengers, C] = ndgrid(Sweep_S, Sweep_W_p, Sweep_v_max, Sweep_passengers, Sweep_C);

W2 = (2.72*S + 0.335 + W_p)*9.81;
W3 = (2.72*S + 0.335 + 0.04*passengers)*9.81;
WL2 = W2./S;
WL3 = W3./S;
T_c = W2.*(0.5*rho*Cd_min*v_max.^2./WL2 + k*2./(rho*v_max.^2).*WL2);
P = T_c.*v_max/(eff_c*eff_e);
v_lof2 = 1.556*sqrt(WL2/(rho*Cl_max));
v_lof3 = 1.556*sqrt(WL3/(rho*Cl_max));
T_lof2 = 1.414*eff_e*P.*polyval(eta_coeff,v_lof2./v_max/1.414)./v_lof2;
T_lof3 = 1.414*eff_e*P.*polyval(eta_coeff,v_lof3./v_max/1.414)./v_lof3;
S_g2 = (W2.*v_lof2.^2)./(2*g*(T_lof2 - 0.5*Cd_min*rho.*S.*v_lof2.^2 - mu*(W2 - 0.5*Cl_max*rho.*S.*v_lof2.^2)));
S_g3 = (W3.*v_lof3.^2)./(2*g*(T_lof3 - 0.5*Cd_min*rho.*S.*v_lof3.^2 - mu*(W3 - 0.5*Cl_max*rho.*S.*v_lof3.^2)));

laps = (v_max.*min(300, C./P)/lap_dist);

M2raw = W_p.*v_max;
M3raw = laps.*passengers./C;
GMraw = (Time2change_config*2 + Time2secure_hatches_doors*2 + ...
    Time2unsecure_hatches_doors*2 + Time2load_m2_payload + ...
    Time2unload_m2_payload + passengers*(Time2loadunload_passenger));
GMScore = GMraw;
M2Score = M2raw;
M3Score = M3raw;

TOFL_constraint = S_g2 > TOFL_max | S_g3 > TOFL_max;
WL_constraint = WL2 > WL_max | WL3 > WL_max;
P_constraint = P > P_max;
lap_constraint = C./P.*v_max < 3.5*lap_dist;
S_constraint = S > S_max;
passenger_constraint = passengers > passenger_max;
all_constraint = TOFL_constraint | P_constraint | lap_constraint | S_constraint | WL_constraint | passenger_constraint;

GMScore(all_constraint) = NaN;
M2Score(all_constraint) = NaN;
M3Score(all_constraint) = NaN;
GMraw = (min(GMScore,[],'all'))./GMraw;
max(M2Score,[],'all')
max(M3Score,[],'all')
M2raw = M2raw/(max(M2Score,[],'all'));
M3raw = M3raw/(max(M3Score,[],'all'));
GMScore = (min(GMScore,[],'all'))./GMScore;
M2Score = M2Score/(max(M2Score,[],'all'));
M3Score = M3Score/(max(M3Score,[],'all'));
Raw = M2raw + M3raw + GMraw;
Score = M2Score + M3Score + GMScore;
[maxScore, idx] = max(Score,[],'all');

[S_i, W_p_i, v_max_i, passenger_i, C_i] = ind2sub(size(Score),idx);
best_config_idx = [S_i, W_p_i, v_max_i, passenger_i, C_i];
best_S = Sweep_S(S_i); 
best_W_p = Sweep_W_p(W_p_i); 
best_v_max = Sweep_v_max(v_max_i);
best_passengers = Sweep_passengers(passenger_i);
best_C = Sweep_C(C_i)/3600;
best_config = [best_S, best_W_p, best_v_max, best_passengers, best_C];
req_Power = P(best_config_idx(1), best_config_idx(2), best_config_idx(3), best_config_idx(4), best_config_idx(5));
req_ST = req_Power*polyval(eta_coeff,0.01)/(0.01*best_config(3))/9.81;
AUW = W2(best_config_idx(1), best_config_idx(2), best_config_idx(3), best_config_idx(4), best_config_idx(5))/9.81;
flight_time = Sweep_C(C_i)/req_Power/60;
req_laps = laps(best_config_idx(1), best_config_idx(2), best_config_idx(3), best_config_idx(4), best_config_idx(5));
v_stall = v_lof2(best_config_idx(1), best_config_idx(2), best_config_idx(3), best_config_idx(4), best_config_idx(5));
T_cruise = T_c(best_config_idx(1), best_config_idx(2), best_config_idx(3), best_config_idx(4), best_config_idx(5));
% Round outputs for excel table
best_S = round(best_S,2);
best_W_p = round(best_W_p,2);
best_v_max = round(best_v_max,2);
best_passengers = round(best_passengers,2);
best_C = round(best_C,2);
req_Power = round(req_Power,2);
req_ST = round(req_ST,2);
AUW = round(AUW,2);
flight_time = round(flight_time,2);
req_laps = round(req_laps,2);
v_stall = round(v_stall,2);
% All of these outputs should go to an excel sheet
T = table(best_S,best_W_p,best_v_max,best_passengers,best_C,req_Power,req_ST,AUW,flight_time,req_laps,v_stall);
T = rows2vars(T); % transpose
filename = 'optimizerV2_outputs.xlsx';
writetable(T,filename)
toc
%% Plotting
c = [0.647 0.078 0.09; 0.424 0.451 0.451];


query = 0.8:0.005:1.2;
v_diff = squeeze(Raw(S_i,W_p_i,:,passenger_i, C_i))/Raw(S_i,W_p_i,v_max_i,passenger_i, C_i);
v_diff = spline(Sweep_v_max(:)/Sweep_v_max(v_max_i), v_diff, query);

P_diff = squeeze(Raw(S_i,W_p_i,:,passenger_i, C_i))/Raw(S_i,W_p_i,v_max_i,passenger_i, C_i);
P_diff = spline(squeeze(P(S_i,W_p_i,:,passenger_i,C_i))/P(S_i,W_p_i,v_max_i,passenger_i, C_i), P_diff, query);

C_diff = squeeze(Raw(S_i,W_p_i,v_max_i,passenger_i,:))/Raw(S_i,W_p_i,v_max_i,passenger_i, C_i);
C_diff = spline(Sweep_C(:)/Sweep_C(C_i), C_diff, query);

figure
plot(query, v_diff)
hold on
plot(Sweep_W_p(:)/Sweep_W_p(W_p_i), squeeze(Raw(S_i,:,v_max_i,passenger_i, C_i))/Raw(S_i,W_p_i,v_max_i,passenger_i, C_i))
plot(Sweep_passengers(:)/Sweep_passengers(passenger_i), squeeze(Raw(S_i,W_p_i,v_max_i,:,C_i))/Raw(S_i,W_p_i,v_max_i,passenger_i, C_i))
plot(query, C_diff)
xlim([0.8 1.2])
hold off

legend('v_{max}', 'Cabinet Weight', 'Passengers', 'Capacity');
xlabel('Fractional Change in Scoring Variable');
ylabel('Fractional Change in Score')

set(0, "DefaultAxesColorOrder", c)

% surf(Sweep_passengers,Sweep_W_p,squeeze(Score(best_config_idx(1),:,best_config_idx(3),:)));