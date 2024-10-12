function init_proplsn_biom
% init function for the standalone propulsion model - populate parameters for
% biom power system - weights x.x kgf

%lipo battery tlu
proplsn.variant_name             = 'biomimetic 740kV, 5 cell lipo, 11x10 prop';
lipo_batt_tlu                    = [3.40 3.65 3.69 3.71 3.73 3.78 3.83 3.90 3.95 4.05 4.20];
proplsn.battery_tlu              = lipo_batt_tlu.*5; %for 5 cell series
proplsn.battery_bkpt_chg_st_dec  = [0.00 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00];
proplsn.V_nom_batt = 5*3.7; %4S lipo
proplsn.max_cap_Ah = 10;
proplsn.internal_resistance = 0.0025;
k_peukert = 1.02; %this is the high discharge loss exponentn for lipo chemistry
proplsn.alpha_peukert = (k_peukert-1)/(2-k_peukert);

% motor params for ???
proplsn.R_sys_ohm = 0.18 + proplsn.internal_resistance + (0.0132*0.5); %Ohm (resistance of the motor plus resistance of the battery plus wiring resistance)
proplsn.K_v = 740;      % rpm/volt
proplsn.K_T_Nm_A = 0.11; % Nm/A <---this was a guess that was iterated on until the efficiency look reasonable
I_0 = 0.5;       % A
proplsn.T_mf_Nm = I_0 * proplsn.K_T_Nm_A;
proplsn.I_xx_rotor_kgm2 = 0.070*(4*.0254)^2; %mass times radius of gyration


% prop tables
% bring in prop data
path = 'UIUCpropDB/';
file = 'grsn_11x10_3002os_4506.txt';
proplsn.prop_dia_in = 11;
proplsn.num_props = 1;
prop_table = importdata([path,file]);
prop_table = prop_table.data;
proplsn.CP_bkpt = prop_table(:,1);
proplsn.CT_tlu  = prop_table(:,2);
proplsn.CP_tlu  = prop_table(:,3);
proplsn.CT_bkpt = proplsn.CP_bkpt;
proplsn.max_J = 1.5;

% write off structure to base
assignin('base','proplsn',proplsn);
end