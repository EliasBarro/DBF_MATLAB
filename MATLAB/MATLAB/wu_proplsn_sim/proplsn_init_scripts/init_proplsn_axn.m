function init_proplsn_axn
% init script for the standalone propulsion model - populate parameters for
% cloudsfly axn - weight is 0.5kg battery

%metadata
proplsn.variant_name     = 'axn 22mm dia motor,nimh 10s2p';
proplsn.num_props = 1;

%nimh battery tlu
num_cells                = 10;
nimh_batt_tlu            = [1.13 1.16 1.18 1.19 1.19 1.20 1.20 1.22 1.24 1.31 1.38];
proplsn.battery_tlu              = nimh_batt_tlu.*num_cells;
proplsn.battery_bkpt_chg_st_dec  = [0.00 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00];
proplsn.V_nom_batt = num_cells*1.2;

%use two elite 1500
proplsn.max_cap_Ah = 3.2;
proplsn.internal_resistance = 0.012*num_cells*0.5;
k_peukert = 1.12; %this is the high discharge loss exponent for nimh chemistry
proplsn.alpha_peukert = (k_peukert-1)/(2-k_peukert);

% set motor params for axn
proplsn.R_sys_ohm = 0.03 + proplsn.internal_resistance + (0.0132*0.5); %Ohm (resistance of the motor plus resistance of the battery plus wiring resistance)
proplsn.K_v = 2200; %rpm/volt
proplsn.K_T_Nm_A = 0.028; %Nm/A <---this was a guess that was iterated on until the efficiency look reasonable
I_0 = 0.4;   %A
proplsn.T_mf_Nm = I_0 * proplsn.K_T_Nm_A;
proplsn.I_xx_rotor_kgm2 = 0.1*(1.1*.0254)^2;

%define path to prop db file
path = 'UIUCpropDB/';
file = 'apcsp_7x6_2936cm_7021.txt';
prop_table = importdata([path,file]);
prop_table = prop_table.data;
proplsn.prop_dia_in = 7; %in

%apply power gain to address mounting considerations
power_gain = 1.0;
prop_table(:,3) = prop_table(:,3).*power_gain;

% output
proplsn.CP_bkpt = prop_table(:,1);
proplsn.CT_tlu  = prop_table(:,2);
proplsn.CP_tlu  = prop_table(:,3);
proplsn.CT_bkpt = proplsn.CP_bkpt;
proplsn.max_J = 1.5;
assignin('base','proplsn',proplsn);

end