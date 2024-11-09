function init_proplsn_revbatwing
% init script for the standalone propulsion model - populate parameters for
% reverse batwing (UNFINISHED)

% get battery data:
current_in = 30;          % estimate 30A pulled
capacity   = 3000;
cells      = 9;
[proplsn.V_nom_batt, proplsn.internal_resistance] = model_nimh_battery(cells, capacity, current_in);

% set motor params for SK3
proplsn.R_sys_ohm = 0.026 + proplsn.internal_resistance + (0.0132*0.5); %Ohm (resistance of the motor plus resistance of the battery plus wiring resistance)
proplsn.K_v = 1500; %rpm/volt
K_t = 0.038; %Nm/A <---this was a guess that was iterated on until the efficiency look reasonable
I_0 = 0.4;   %A
torque_mf = I_0 * K_t;

%define path to prop db file
path = 'UIUCpropDB/';
file = 'grcsp_9x5_kt0963_6922.txt';
prop_diameter = 9; %in
prop_filename = [path, file];
proplsn.num_props = 1;

%apply power gain to address mounting considerations
power_gain = 1.5;
prop_table(:,3) = prop_table(:,3).*power_gain;

% outputs
proplsn.CP_bkpt = prop_table(:,1);
proplsn.CT_tlu  = prop_table(:,2);
proplsn.CP_tlu  = prop_table(:,3);
proplsn.CT_bkpt = proplsn.CP_bkpt;
proplsn.max_J = 1.5;
assignin('base','proplsn',proplsn);

end