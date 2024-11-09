%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first_order_motor_model_static
% Author: Austin Stover
% Date: October-November 2018
% Plot data on a propulsion configuration, based on prop data from the UIUC
% database and Mark Drela's first order DC motor model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%TODO: Determine electrical efficiency at fixed thrust
% - Find required voltage for each Kv and desired T and run sim on required V for each Kv

T_des_kgf_list = 0.2:1:3.2; %4:1:8; %Desired thrust

% Aero Params
V = 0; % 25; %Estimated velocity (m/s). 0 requires static thrust data
rho = 1.16; %Air density (kg/m^3)

%Propulsion Params
num_series = 6; %23;
num_parallel = 1; %2;
K_v_list = 50:10:850; %Motor Kv (RPM/V)
I_0 = 1; %Motor no-load current (A)
R_m = 0.1; %Motor internal resistance (Ohms)

v_nom_cell = 3.4; %1.2; %NiMh
R_cell = 0.004; 
R_wires = 0.0047*0.5;
I_cell_max_discharge = 30; %A

I_max = I_cell_max_discharge*num_parallel;
v_nom = v_nom_cell*num_series; %Voltage
R_circuit = R_cell*num_series/num_parallel + R_wires;
T_des = T_des_kgf_list*9.80665;


path = 'UIUCpropDB/';
file = 'rusp_11x4_static_2952os.txt'; %'15x10_static_QProp.txt'; %'15x10at8700rpmQProp.txt'; %'rusp_11x4_static_2952os.txt';
prop_table = importdata([path,file]);
prop_table = prop_table.data;
prop_dia_in = 15; %Use the actual prop dia.
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

mult_num_pts = 50; %Multiply the number of datapoints by this constant w/ spline interpolation
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

I_list = zeros(length(RPM),1,length(K_v_list)); %TODO: Add velocity of plane as y in array
v_list = zeros(length(RPM),1,length(K_v_list));
Q_m_list = zeros(length(RPM),1,length(K_v_list));
eta_m_list = zeros(length(RPM),1,length(K_v_list));
RPMIntersec_list = zeros(1,length(K_v_list));
T_intersec_list = zeros(1,length(K_v_list));
I_intersec_list = zeros(1,length(K_v_list));
v_nom_des_list = zeros(1,length(K_v_list),length(T_des_kgf_list));
for i = 1:length(K_v_list)
	K_v = K_v_list(i);
	K_v_rad = 2*pi/60*K_v; %Kv in rad/s/V
	Omega = 2*pi*n; %rad/s

	% Current
	I = (v_nom - Omega./K_v_rad)./(R_m + R_circuit);

	% Voltage across motor v = v_nom - R_circuit*i
	v = v_nom - R_circuit.*I;

	Q_m = ((v - Omega./K_v_rad)./R_m - I_0)./K_v_rad; %Motor torque

	eta_m = (1 - I_0*R_m./(v - Omega./K_v_rad)) .* Omega./(v.*K_v_rad); %Motor efficiency

	% Find intersec of lines connecting pts btwn QIntersecInd, QIntersecInd+1
	QIntersecInd = find(diff(Q - Q_m > 0),1); % Q > Q_m; % Find the 2 points spanning the intersection of Q and Q_m
	if ~isempty(QIntersecInd) %If there is an intersection after the first val of RPM
		Q_slope = (Q(QIntersecInd) - Q(QIntersecInd+1))/(RPM(QIntersecInd) - RPM(QIntersecInd+1));
		Q_m_slope = (Q_m(QIntersecInd) - Q_m(QIntersecInd+1))/(RPM(QIntersecInd) - RPM(QIntersecInd+1));
		RPMIntersec = ((Q_slope - Q_m_slope)*RPM(QIntersecInd) - Q(QIntersecInd) + Q_m(QIntersecInd))...
		/(Q_slope - Q_m_slope);
	T_intersec_list(i) = T(QIntersecInd);
	I_intersec_list(i) = I(QIntersecInd);
	else
		RPMIntersec = RPM(end);
	end
	
	for j = 1:length(T_des_kgf_list)
		[~,T_des_ind] = min(abs(T - T_des(j))); %Find closest RPM pt to intersection
		v_req = (K_v_rad*Q(T_des_ind) + I_0)*R_m + Omega(T_des_ind)/K_v_rad; %Voltage required for the desired thrust
		I_req = (v_req - Omega(T_des_ind)./K_v_rad)./(R_m + R_circuit);
		v_nom_des = v_req + (R_circuit/R_m)*(v_req - Omega(T_des_ind)/K_v_rad);
% 		bat_maxed = I_req > I_max;
		v_nom_des_list(1,i,j) = v_nom_des;%*(1-bat_maxed);
	end
	
	I_list(:,1,i) = I;
	v_list(:,1,i) = v;
	Q_m_list(:,1,i) = Q_m;
	eta_m_list(:,1,i) = eta_m;
	RPMIntersec_list(1,i) = RPMIntersec;
end

% fprintf('Recommended voltage for desired %gkgf thrust: %g\n',T_des_kgf,v_nom_des);
% fprintf('Recommended battery pack size: %gS\n',ceil(v_nom_des/v_nom_cell));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure; %Plot CP and CT vs RPM
% hold on;
% yyaxis left;
% plot(RPM, CT);
% xlabel('RPM');
% ylabel('C_T');
% yyaxis right;
% plot(RPM,CP);
% ylabel('C_P');
% hold off;

figure;
hold on;
for j = 1:length(T_des_kgf_list)
	plot(K_v_list,v_nom_des_list(1,:,j)/v_nom_cell);
end
legend("T = " + string(T_des_kgf_list) + " Kgf");
ylabel('Series Cell Count');
xlabel('K_v');
if(V ~= 0)
	title(['V = ' num2str(V)])
else
	title('Static')
end

figure;
title(['V = ' num2str(V)])
yyaxis left;
plot(K_v_list, T_intersec_list./9.80665);
ylabel('Dynamic Thrust (kgf)');
xlabel('K_v');
yLimits = ylim;
ylim([0,max(realmin,yLimits(2))]);
yyaxis right;
plot(K_v_list, (1000*T_intersec_list./9.80665) ./ (I_intersec_list.*v_nom));
ylabel('Efficiency (gf/W)');
xlabel('K_v');
yLimits = ylim;
ylim([0,max(realmin,yLimits(2))]);


%Plot x-limits:
RPMLim = [min(RPM),max(RPM)]; %[min(RPM),max(RPM)];

figure;
if(V ~= 0)
	ipanel(@(K_v_ind)motor_data_plot(V, K_v_list(floor(K_v_ind)),RPM,RPMLim, ...
	RPMIntersec_list(1,floor(K_v_ind)),T,Q,Q_m_list(:,1,floor(K_v_ind)), ...
		I_list(:,1,floor(K_v_ind)),v_list(:,1,floor(K_v_ind)),I_max, ...
		eta_m_list(:,1,floor(K_v_ind)),eta), ...
   {'slider','K_v Index',{1,length(K_v_list)+0.5,1}},'MinControlWidth',300,'LabelWidth',100);
else %For static thrust sims, eta does not exist
	ipanel(@(K_v_ind)motor_data_plot(V, K_v_list(floor(K_v_ind)),RPM,RPMLim, ...
			RPMIntersec_list(1,floor(K_v_ind)),T,Q,Q_m_list(:,1,floor(K_v_ind)), ...
				I_list(:,1,floor(K_v_ind)),v_list(:,1,floor(K_v_ind)),I_max, ...
				eta_m_list(:,1,floor(K_v_ind))), ...
		   {'slider','K_v Index',{1,length(K_v_list)+0.5,1}},'MinControlWidth',300,'LabelWidth',100);
end
	