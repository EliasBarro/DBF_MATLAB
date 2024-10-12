%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% motor_model_static_and_dynamic
% Author: Austin Stover
% Date: October-November 2018
% Plot data on a propulsion configuration, based on prop data from the UIUC
% database and Mark Drela's first order DC motor model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

stat_T_des_kgf = 40; %Desired static thrust
% Aero Params
V = 0; %Estimated velocity (m/s) for dynamic thrust
rho = 1.16; %Air density (kg/m^3)

%Propulsion Params
K_v_list = 5:10:850; %Motor Kv (RPM/V)
I_0 = 1; %Motor no-load current (A)
R_m = 0.1; %Motor internal resistance (Ohms)
v_nom_cell = 3.4; %1.2;


stat_T_des = stat_T_des_kgf*9.80665;

path = 'UIUCpropDB/';
%Dynamic prop data
dyn_file = 'apce_19x12_jb1084_3007.txt'; %'15x10at8700rpmQProp.txt'; %'rusp_11x4_static_2952os.txt';
dyn_prop_table = importdata([path,dyn_file]);
dyn_prop_table = dyn_prop_table.data;
%Static prop data
stat_file = 'apce_19x12_static_jb1078.txt'; %15x10_static_QProp.txt'; %'rusp_11x4_static_2952os.txt';
stat_prop_table = importdata([path,stat_file]);
stat_prop_table = stat_prop_table.data;

prop_dia_in = 15; %Prop Diameter
n_blades = 3;
mult_num_pts = 50; %Multiply the number of datapoints by this constant w/ spline interpolation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prop Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prop tables
% bring in prop data
% apply power gain
power_gain = n_blades/2;
dyn_prop_table(:,3) = dyn_prop_table(:,3).*power_gain;
stat_prop_table(:,3) = stat_prop_table(:,3).*power_gain;

% apply thrust gain
thrust_gain = max(1,0.9*n_blades/2);
dyn_prop_table(:,2) = dyn_prop_table(:,2).*thrust_gain;
stat_prop_table(:,2) = stat_prop_table(:,2).*thrust_gain;

dyn_sz = size(dyn_prop_table,1);
stat_sz = size(stat_prop_table,1);

%Dynamic data

% stretch J for more pitch
J_gain = 1; % No modifications to pitch
dyn_prop_table(:,1) = dyn_prop_table(:,1).*J_gain;

dyn_eta = dyn_prop_table(:,4); % Aerodynamic efficiency
dyn_CT  = dyn_prop_table(:,2); % Coefficient of thrust
dyn_CP  = smooth(dyn_prop_table(:,3)); % Coefficient of power
dyn_J = dyn_prop_table(:,1); % Advance Ratio J = V/(nD); Rotations per sec n = RPM/60

dyn_eta = interp1(1:dyn_sz,dyn_eta,1:1/mult_num_pts:dyn_sz,'spline');
dyn_CT = interp1(1:dyn_sz,dyn_CT,1:1/mult_num_pts:dyn_sz,'spline');
dyn_CP = interp1(1:dyn_sz,dyn_CP,1:1/mult_num_pts:dyn_sz,'spline');
dyn_J = interp1(1:dyn_sz,dyn_J,1:1/mult_num_pts:dyn_sz,'spline');

D = prop_dia_in*0.0254; %Dia. (m)
dyn_RPM = 60*V./(dyn_J.*D);

%Static Data
stat_RPM = stat_prop_table(:,1);
stat_CT  = smooth(stat_prop_table(:,2)); % Coefficient of thrust
stat_CP  = smooth(stat_prop_table(:,3)); % Coefficient of power

stat_RPM = interp1(1:stat_sz,stat_RPM,1:1/mult_num_pts:stat_sz,'spline');
stat_CT = interp1(1:stat_sz,stat_CT,1:1/mult_num_pts:stat_sz,'spline');
stat_CP = interp1(1:stat_sz,stat_CP,1:1/mult_num_pts:stat_sz,'spline');


dyn_n = dyn_RPM./60; % n is rot/s
stat_n = stat_RPM./60;

% Coefficient of thrust CT = T/(rho n^2 D^4)  See McCormick Aerodynamics, Aeronautics, and Flight Mechanics for more info
dyn_T = dyn_CT.*(rho .* dyn_n.^2 .* D^4);
stat_T = stat_CT.*(rho .* stat_n.^2 .* D^4);

% Coefficient of power CP = P/(rho n^3 D^5)
dyn_P = dyn_CP.*(rho .* dyn_n.^3 .* D^5);
stat_P = stat_CP.*(rho .* stat_n.^3 .* D^5);

% P = omega*Q, where Q is torque and omega is rad/s
dyn_Q = dyn_P./(2*pi*dyn_n); %Prop torque
stat_Q = stat_P./(2*pi*stat_n); %Prop torque



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Motor Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (See Drela's motor-prop matching notes, http://web.mit.edu/drela/Public/web/qprop/motorprop.pdf
% Current through motor: i(Omega,v) = (v - Omega/Kv)/R_m
% Motor RPM: Omega(i,v) = (v - iR_m)Kv
% Motor torque: Qm(Omega,v) = [(v-Omega/Kv)/R_m - i0]/Kv
% Motor efficiency: eta_m(i,v) = P_shaft/P_elec = (1-i0/i)(1 - iR_m/v)
% Assuming K_v is in rad/s/V, Omega is in rad/s


%Static model (at desired static thrust)
for i = 1:length(K_v_list)
	K_v = K_v_list(i);
	K_v_rad = 2*pi/60*K_v; %Kv in rad/s/V
	stat_Omega = 2*pi*stat_n; %rad/s
	[~,stat_T_des_ind] = min(abs(stat_T - stat_T_des)); %Find RPM pt where T is closest to T_des
	stat_v_req = (K_v_rad*stat_Q(stat_T_des_ind) + I_0)*R_m + stat_Omega(stat_T_des_ind)/K_v_rad; %Voltage required for the desired thrust
	stat_v = stat_v_req;
	stat_Q_m = ((stat_v - stat_Omega./K_v_rad)./R_m - I_0)./K_v_rad; %Motor torque
	stat_eta_m = (1 - I_0*R_m./(stat_v - stat_Omega./K_v_rad)) .* stat_Omega./(stat_v.*K_v_rad); %Motor efficiency
	
	% Find ss values
	% Find intersec of lines connecting pts btwn QIntersecInd, QIntersecInd+1
	[~,QIntersecInd] = min(abs(stat_Q - stat_Q_m)); %Find the pt of intersection btwn Q, Q_m
	
	stat_P_mech = stat_Omega(QIntersecInd)*stat_Q_m(QIntersecInd);
	stat_P_elec = stat_P_mech/stat_eta_m(QIntersecInd);
	
	stat_v_req_list(i) = stat_v_req;
	stat_T_list(i) = stat_T(QIntersecInd);
	stat_P_elec_list(i) = stat_P_elec;
	
	%I_list(i) = I(QIntersecInd);
	%v_nom_list(i) = v_nom_req;
	%eta_m_list(i) = eta_m(QIntersecInd);
	%v_m_list(i) = v;
	%I_list(i) = I;
end

%Dynamic model (at required voltage for static thrust)
for i = 1:length(K_v_list)
	K_v = K_v_list(i);
	K_v_rad = 2*pi/60*K_v; %Kv in rad/s/V
	dyn_Omega = 2*pi*dyn_n; %rad/s
	dyn_v = stat_v_req_list(i); %Required voltage for static thrust at this K_v
	dyn_Q_m = ((dyn_v - dyn_Omega./K_v_rad)./R_m - I_0)./K_v_rad; %Motor torque
	dyn_eta_m = (1 - I_0*R_m./(dyn_v - dyn_Omega./K_v_rad)) .* dyn_Omega./(dyn_v.*K_v_rad); %Motor efficiency
	
	% Find ss values
	% Find intersec of lines connecting pts btwn QIntersecInd, QIntersecInd+1
	[~,QIntersecInd] = min(abs(dyn_Q - dyn_Q_m)); %Find the pt of intersection btwn Q, Q_m
	
	dyn_P_mech = dyn_Omega(QIntersecInd)*dyn_Q_m(QIntersecInd);
	dyn_P_elec = dyn_P_mech/dyn_eta_m(QIntersecInd);
	
	dyn_T_list(i) = dyn_T(QIntersecInd);
	dyn_P_elec_list(i) = dyn_P_elec;
	
end

stat_elec_eff_gf_W = (stat_T_list/9.8055*1000)./stat_P_elec_list;
dyn_elec_eff_gf_W = (dyn_T_list/9.8055*1000)./dyn_P_elec_list;

figure;
hold on;
plot(K_v_list, stat_elec_eff_gf_W);
plot(K_v_list, dyn_elec_eff_gf_W);
title(sprintf('Static Thrust = %g kgf',stat_T_des_kgf));
xlabel('K_v');
ylabel('Electrical Efficiency [gf/W]');
yLimits = ylim;
ylim([0,max(realmin,yLimits(2))]);
%ylim([0,7]);
xlim([min(K_v_list), max(K_v_list)]);
legend('Static',sprintf('V = %g m/s',V));

% figure;
% plot(K_v_list, dyn_T_list/9.8055);
% title(sprintf('Dynamic Thrust for Static Thrust = %g kgf',stat_T_des_kgf));
% xlabel('K_v');
% ylabel('Thrust kgf');
% yLimits = ylim;
% ylim([0,max(realmin,yLimits(2))]);

figure;
plot(K_v_list, stat_v_req_list/v_nom_cell);
title('Series Cell Count');
ylim([0,60]);