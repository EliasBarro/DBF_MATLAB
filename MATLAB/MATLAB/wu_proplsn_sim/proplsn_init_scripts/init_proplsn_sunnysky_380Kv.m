function init_proplsn_sunnysky_380Kv
% init function for the standalone propulsion model - populate parameters for
% Motor: https://sunnyskyusa.com/collections/x-motors/products/sunnysky-x4130-brushless-motors?variant=45674524431
% Max. continuous power:	2350W
% Max. continuous current:	79A
% Mass:						416g

% NiMH battery tlu
num_series                = 24; %Turnigy 5000mAh: 0.070 kg/cell
num_parallel			  = 2;
nimh_batt_tlu            = [1.13 1.16 1.18 1.19 1.19 1.20 1.20 1.22 1.24 1.31 1.38];
proplsn.battery_tlu              = nimh_batt_tlu.*num_series;
proplsn.battery_bkpt_chg_st_dec  = [0.00 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00];
proplsn.V_nom_batt = num_series*1.2;

%Use Turnigy 5000s in series. Turnigy 5000s: https://hobbyking.com/en_us/turnigy-sub-c-1-2v-5000mah-high-power-series-nimh-single-cell.html
cell_max_cap_Ah = 5;
cell_internal_resistance = 0.008; %Found based on actual tested motor thrust with NiMH
proplsn.max_cap_Ah = cell_max_cap_Ah*num_parallel;
proplsn.internal_resistance = cell_internal_resistance*num_series/num_parallel;
k_peukert = 1.12; %this is the high discharge loss exponent for nimh chemistry
proplsn.alpha_peukert = (k_peukert-1)/(2-k_peukert);

% %lipo battery tlu
% num_series = 12;
% num_parallel = 2;
% lipo_batt_tlu            = [ 3.40 3.65 3.69 3.71 3.73 3.78 3.83 3.90 3.95 4.05 4.20];
% proplsn.battery_tlu              = lipo_batt_tlu.*num_series;
% proplsn.battery_bkpt_chg_st_dec  = [0.00 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00];
% proplsn.V_nom_batt = num_series*3.7; % LiPo nominal voltage = 3.7V
% cell_max_cap_Ah = 5;
% cell_internal_resistance = 0.002;
% 
% proplsn.max_cap_Ah = cell_max_cap_Ah*num_parallel;
% proplsn.internal_resistance = cell_internal_resistance*num_series/num_parallel;
% k_peukert = 1.02;           % this is the high discharge loss exponentn for lipo chemistry
% proplsn.alpha_peukert = (k_peukert-1)/(2-k_peukert);

path = 'UIUCpropDB/';
file = 'apce_17x12_rd1098_3407.txt';% '15x10at8700rpmQProp.txt'; %Propeller: https://www.apcprop.com/product/15x8e/
prop_table = importdata([path,file]);
prop_table = prop_table.data;
proplsn.prop_dia_in = 18; %Scales prop: can tweak this by small amounts (1" for >12" dia.; 0.5" for <12" dia. props)
propeller_mass = 0.0439;
n_blades = 2;
proplsn.num_props = 2;

% motor params for Sunnysky 380Kv
proplsn.R_sys_ohm = 0.021 + proplsn.internal_resistance + (0.0047*0.5); %Ohm (resistance of the motor plus resistance of the battery plus wiring resistance)
proplsn.K_v = 323.9;      % Kv, rpm/volt
proplsn.K_T_Nm_A = 1/proplsn.K_v * 60 / (2*pi);
I_0 = 2.1;              % No-load current A   % ^Basically, make sure that mechanical power is somewhere around 0.9 times electrical power
proplsn.T_mf_Nm = I_0 * proplsn.K_T_Nm_A;
proplsn.I_xx_rotor_kgm2 = propeller_mass*(1/12)*(proplsn.prop_dia_in/2*0.0254)^2; %mass times radius of gyration
stability_factor        = 4; %Higher -> less wiggles in RPM in proplsn_table
proplsn.I_xx_rotor_kgm2 = proplsn.I_xx_rotor_kgm2*stability_factor;


% apply power gain
power_gain = n_blades/2;
prop_table(:,3) = prop_table(:,3).*power_gain;

% apply thrust gain
thrust_gain = max(1,0.9*n_blades/2);
prop_table(:,2) = prop_table(:,2).*thrust_gain;

% stretch J for more pitch
J_gain = 12/12;
prop_table(:,1) = prop_table(:,1).*J_gain;


% proplsn table options. Run `get_proplsn_table(false, true)` to determine 
% these vals. Good vals allow the thrust, current, and rpm to settle out
% after each step.
%settle_times define the propulsion table timestep lengths in
%proplsn_table. They must be increasing from left to right.
proplsn.settle_times = [0.25, 1.5, 1.5]; 
proplsn.bp_dens = linspace(0.9,1.225,6);
proplsn.bp_u_fs = linspace(0,50,25);
proplsn.bp_thr  = linspace(0,1,6);
proplsn.dt      = 0.001; %The simulation timestep. Should be smaller than 0.01 (100 hertz) if oscillations

% output
proplsn.CP_bkpt = prop_table(:,1);
proplsn.CT_tlu  = prop_table(:,2); % Coefficient of thrust
proplsn.CP_tlu  = smooth(prop_table(:,3)); % Coefficient of power
proplsn.CT_bkpt = proplsn.CP_bkpt;
proplsn.max_J = max(prop_table(:,1));
assignin('base','proplsn',proplsn);

end