function init_proplsn_v4_Flight1a
% Second modification of v2: With a reasonable estimate for battery and
% less dynamic thrust
% capacity and mass

% Motor: ?
% Prop: ?
% Mass: ~1.5kg (2kgf max)

% lipo battery tlu
num_series = 12;
num_parallel = 1;
lipo_batt_tlu            = [ 3.40 3.65 3.69 3.71 3.73 3.78 3.83 3.90 3.95 4.05 4.20];
proplsn.battery_tlu              = lipo_batt_tlu.*num_series;
proplsn.battery_bkpt_chg_st_dec  = [0.00 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00];
proplsn.V_nom_batt = num_series*3.7; % LiPo nominal voltage = 3.7V
cell_max_cap_Ah = 5;
cell_internal_resistance = 0.001125;

proplsn.max_cap_Ah = cell_max_cap_Ah*num_parallel;
proplsn.internal_resistance = cell_internal_resistance*num_series/num_parallel;
k_peukert = 1.02;           % this is the high discharge loss exponent for lipo chemistry
proplsn.alpha_peukert = (k_peukert-1)/(2-k_peukert);

% prop tables
% bring in prop data
path = 'UIUCpropDB/';
file = 'apce_11x8_pg0524_6000.txt'; %'apce_9x6_rd0994_6715.txt';%'apce_9x7.5_2551rd_6014.txt';%'apce_8x8_2812rd_6420.txt'; %'apce_8x6_2804rd_6712.txt';
prop_table = importdata([path,file]);
prop_table = prop_table.data;
proplsn.prop_dia_in = 12;
propeller_mass = 0.040;
n_blades = 2;
proplsn.num_props = 2;

% motor params for EF1300 Probably this: https://hobbyking.com/en_us/ntm-prop-drive-series-ef-1-pylon-racing-motor-1300kv-930w-v2.html
proplsn.R_sys_ohm = 0.02 + proplsn.internal_resistance + (0.0132*0.5); %Ohm (resistance of the motor plus resistance of the battery plus wiring resistance)
proplsn.K_v = 323.9;      % rpm/volt
proplsn.K_T_Nm_A = (1/proplsn.K_v)*(60/(2*pi));
I_0 = 2.1;               % A
proplsn.T_mf_Nm = I_0 * proplsn.K_T_Nm_A;
proplsn.I_xx_rotor_kgm2 = propeller_mass*(1.75*.0254)^2; %mass times radius of gyration
stability_factor        = 1.3;
proplsn.I_xx_rotor_kgm2 = proplsn.I_xx_rotor_kgm2*stability_factor;


%apply power gain to address mounting considerations
power_gain_n_blades = n_blades/2;
power_gain = power_gain_n_blades*1.02; %slight reduction, since fairly slim, but there is pusher/wing interactions
prop_table(:,3) = prop_table(:,3).*power_gain;

% apply thrust gain
thrust_gain = max(1,0.9*n_blades/2);
prop_table(:,2) = prop_table(:,2).*thrust_gain;

% stretch J for more pitch
J_gain = 8/8;
prop_table(:,1) = prop_table(:,1).*J_gain;

% proplsn table options
proplsn.settle_times = [0.15, 2, 2];
proplsn.bp_dens = linspace(0.9,1.225,6);
proplsn.bp_u_fs = linspace(0,50,25);
proplsn.bp_thr  = linspace(0,1,6);
proplsn.dt      = 0.01;

% output
proplsn.CP_bkpt = prop_table(:,1);
proplsn.CT_tlu  = prop_table(:,2);
proplsn.CP_tlu  = prop_table(:,3);
proplsn.CT_bkpt = proplsn.CP_bkpt;
proplsn.max_J = max(prop_table(:,1));
assignin('base','proplsn',proplsn);

end