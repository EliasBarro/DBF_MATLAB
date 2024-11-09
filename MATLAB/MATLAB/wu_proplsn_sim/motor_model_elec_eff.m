%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first_order_motor_elec_eff
% Author: Austin Stover
% Date: October-November 2018
% Plot data on a propulsion configuration, based on prop data from the UIUC
% database or QProp and Mark Drela's first order DC motor model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

T_des_kgf = 0.66; %Desired thrust
% Aero Params
V = 0; %Estimated velocity (m/s) for dynamic thrust. V = 0 requires static thrust data
rho = 1.16; %Air density (kg/m^3)

%Propulsion Params
K_v_list = 0:10:850; %Motor Kv (RPM/V)
I_0 = 1; %Motor no-load current (A)
R_m = 0.1; %Motor internal resistance (Ohms)
T_des = T_des_kgf*9.80665;

path = 'UIUCpropDB/';
file = 'rusp_11x4_static_2952os.txt'; %'15x10_static_QProp.txt';
mult_num_pts = 50; %Multiply the number of datapoints by this constant w/ spline interpolation
prop_table = importdata([path,file]);
prop_table = prop_table.data;
prop_dia_in = 11; %Use the actual prop dia.
n_blades = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prop Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prop tables
% bring in prop data
% apply power gain
power_gain = n_blades/2;
prop_table(:,3) = prop_table(:,3).*power_gain;

% apply thrust gain
thrust_gain = max(1,0.9*n_blades/2);
prop_table(:,2) = prop_table(:,2).*thrust_gain;

sz = size(prop_table,1);

if(V ~= 0)
	% stretch J for more pitch
	J_gain = 1; % No modifications to pitch
	prop_table(:,1) = prop_table(:,1).*J_gain;

	eta = prop_table(:,4); % Aerodynamic efficiency
	CT  = prop_table(:,2); % Coefficient of thrust
	CP  = smooth(prop_table(:,3)); % Coefficient of power
	J = prop_table(:,1); % Advance Ratio J = V/(nD); Rotations per sec n = RPM/60
	
	eta = interp1(1:sz,eta,1:1/mult_num_pts:sz,'spline');
	CT = interp1(1:sz,CT,1:1/mult_num_pts:sz,'spline');
	CP = interp1(1:sz,CP,1:1/mult_num_pts:sz,'spline');
	J = interp1(1:sz,J,1:1/mult_num_pts:sz,'spline');

	D = prop_dia_in*0.0254; %Dia. (m)
	RPM = 60*V./(J.*D);
	
else %Static thrust UIUC data
	RPM = prop_table(:,1);
	CT  = smooth(prop_table(:,2)); % Coefficient of thrust
	CP  = smooth(prop_table(:,3)); % Coefficient of power
	
	RPM = interp1(1:sz,RPM,1:1/mult_num_pts:sz,'spline');
	CT = interp1(1:sz,CT,1:1/mult_num_pts:sz,'spline');
	CP = interp1(1:sz,CP,1:1/mult_num_pts:sz,'spline');
	
	D = prop_dia_in*0.0254; %Dia. (m)
end

n = RPM./60; % n is rot/s

% Coefficient of thrust CT = T/(rho n^2 D^4)  See McCormick Aerodynamics, Aeronautics, and Flight Mechanics for more info
T = CT.*(rho .* n.^2 .* D^4);

% Coefficient of power CP = P/(rho n^3 D^5)
P = CP.*(rho .* n.^3 .* D^5);

% P = omega*Q, where Q is torque and omega is rad/s
Q = P./(2*pi*n); %Prop torque



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Motor Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (See Drela's motor-prop matching notes, http://web.mit.edu/drela/Public/web/qprop/motorprop.pdf
% Current through motor: i(Omega,v) = (v - Omega/Kv)/R_m
% Motor RPM: Omega(i,v) = (v - iR_m)Kv
% Motor torque: Qm(Omega,v) = [(v-Omega/Kv)/R_m - i0]/Kv
% Motor efficiency: eta_m(i,v) = P_shaft/P_elec = (1-i0/i)(1 - iR_m/v)
% Assuming K_v is in rad/s/V, Omega is in rad/s

for i = 1:length(K_v_list)
	K_v = K_v_list(i);
	K_v_rad = 2*pi/60*K_v; %Kv in rad/s/V
	Omega = 2*pi*n; %rad/s
	[~,T_des_ind] = min(abs(T - T_des)); %Find RPM pt where T is closest to T_des
	v_req = (K_v_rad*Q(T_des_ind) + I_0)*R_m + Omega(T_des_ind)/K_v_rad; %Voltage required for the desired thrust
	v = v_req;
	Q_m = ((v - Omega./K_v_rad)./R_m - I_0)./K_v_rad; %Motor torque
	eta_m = (1 - I_0*R_m./(v - Omega./K_v_rad)) .* Omega./(v.*K_v_rad); %Motor efficiency
	
	% Find ss values
	% Find intersec of lines connecting pts btwn QIntersecInd, QIntersecInd+1
	[~,QIntersecInd] = min(abs(Q - Q_m)); %Find the pt of intersection btwn Q, Q_m
	
	P_mech = Omega(QIntersecInd)*Q_m(QIntersecInd);
	P_elec = P_mech/eta_m(QIntersecInd);
	
	T_list(i) = T(QIntersecInd);
	P_elec_list(i) = P_elec;
end

elec_eff_gf_W = (T_list/9.8055*1000)./P_elec_list;

figure;
plot(K_v_list, elec_eff_gf_W);
xlabel('K_v')
ylabel('Electrical Efficiency [gf/W]')
yLimits = ylim;
ylim([0,max(realmin,yLimits(2))]);