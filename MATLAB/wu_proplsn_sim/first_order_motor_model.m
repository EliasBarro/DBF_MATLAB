%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first_order_motor_model
% Author: Austin Stover
% Date: October 2018
% Plot data on a propulsion configuration, based on prop data from the UIUC
% database and Mark Drela's first order DC motor model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T_des_kgf = 5; %Desired thrust

% Aero Params
V = 23; %Estimated velocity (m/s)
rho = 1.16; %Air density (kg/m^3)

%Propulsion Params
num_series = 22;
num_parallel = 2;
K_v = 380; %Motor Kv (RPM/V)
I_0 = 3.69; %Motor no-load current (A)
R_m = 0.009; %Motor internal resistance (Ohms)

v_nom_cell = 1.2; %NiMh
R_cell = 0.002;
R_wires = 0.0047*0.5;

v_nom = v_nom_cell*num_series; %Voltage
R_circuit = R_cell*num_series/num_parallel + R_wires;
T_des = T_des_kgf*9.80665;


path = 'UIUCpropDB/';
file = '15x10at8700rpmQProp.txt'; %Propeller: https://www.apcprop.com/product/14x12/
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

% stretch J for more pitch
J_gain = 8/10; % Pitch at 8
prop_table(:,1) = prop_table(:,1).*J_gain;

eta = prop_table(:,4); % Aerodynamic efficiency
CT  = prop_table(:,2); % Coefficient of thrust
CP  = smooth(prop_table(:,3)); % Coefficient of power
J = prop_table(:,1); % Advance Ratio J = V/(nD); Rotations per sec n = RPM/60

D = prop_dia_in*0.0254; %Dia. (m)
RPM = 60*V./(J.*D);
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
else
	RPMIntersec = RPM(end);
end

[~,T_des_ind] = min(abs(T - T_des)); %Find closest RPM pt to intersection
v_req = (K_v_rad*Q(T_des_ind) + I_0)*R_m + Omega(T_des_ind)/K_v_rad; %Voltage required for the desired thrust
v_nom_des = v_req + (R_circuit/R_m)*(v_req - Omega(T_des_ind)/K_v_rad);
fprintf('Recommended voltage for desired %gkgf thrust: %g\n',T_des_kgf,v_nom_des);
fprintf('Recommended battery pack size: %gS\n',ceil(v_nom_des/v_nom_cell));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RPMLim = [min(RPM),max(RPM)];

% figure;
% subplot(3,1,1);
% plot(J, eta);
% xlabel('J');
% ylabel('\eta');

NUM_SUBPLOTS_Y = 4;
NUM_SUBPLOTS_X = 1;

figure;
subplot(NUM_SUBPLOTS_Y,NUM_SUBPLOTS_X,1);
hold on;
eta_plot = plot(RPM, eta);
eta_m_plot = plot(RPM, eta_m);
hold off;
xlim(RPMLim);
ylim([0,1]);
line([RPMIntersec, RPMIntersec], ylim,'Color','black','LineStyle','--')
%xlabel('RPM');
ylabel('\eta');
legend([eta_plot,eta_m_plot], {'\eta_{prop}', '\eta_{motor}'},'location','best');

subplot(NUM_SUBPLOTS_Y,NUM_SUBPLOTS_X,2);
plot(RPM, T/9.80665);
xlim(RPMLim);
line(xlim,[T_des_kgf,T_des_kgf],'Color','black','LineStyle','--');
line([RPMIntersec, RPMIntersec], ylim,'Color','black','LineStyle','--')
%xlabel('RPM');
ylabel('T (kgf)');

subplot(NUM_SUBPLOTS_Y,NUM_SUBPLOTS_X,3);
hold on;
prop_Q_plot = plot(RPM, Q);
xlim(RPMLim);
%xlabel('RPM');
ylabel('Q (Nm)');

subplot(NUM_SUBPLOTS_Y,NUM_SUBPLOTS_X,3);
motor_Q_plot = plot(RPM, Q_m);
xlim(RPMLim)
Q_plot_y_lim = ylim;
ylim([0,Q_plot_y_lim(2)])
line([RPMIntersec, RPMIntersec], ylim,'Color','black','LineStyle','--')
legend([prop_Q_plot, motor_Q_plot], {'Q', 'Q_m'},'location','best');
hold off;

subplot(NUM_SUBPLOTS_Y,NUM_SUBPLOTS_X,4);
yyaxis left;
plot(RPM, I);
xlim(RPMLim);
ylim([0,1.5*max(I)])
line([RPMIntersec, RPMIntersec], ylim,'Color','black','LineStyle','--')
xlabel('RPM');
ylabel('I (A)');
yyaxis right;
plot(RPM,I.*v);
ylabel('P_{elec} (W)');
ylim([0,1.5*max(I.*v)])
