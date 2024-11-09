max_power = 1000;
S = 0.5;
AUW = 7;
Cd_min = 0.086;

Cd_i = 0.015; %Induced drage coefficient (taken from GTech 2016 plane)
Cl_to = 0.75; %lift at 0 AoA (taken from GTech 2016 plane)
Cl_max = 1.44; %Max lift at optimal AoA
p = 1.225; %Air density
mu = 0.05; %Ground roll friction coefficient
e_eff = 0.9; %Electrical Efficiency (P_mech/P_elec)
p_eff = 0.83; %Propulsive Efficiency (P_pro/P_mech)
p_eff_to = 0.35; %Propulsive Efficiency at takeoff velocities (P_pro/P_mech) -> propulsive efficeincy decreases with foward velocity (See Snorri pg. 619)
corr = 1.5; %Correction factor to make equation consistent with 2021

v_lof = 1.556*sqrt(AUW*9.81/(p*S*Cl_max));
L = 0.5*p*S*Cl_to*(v_lof^2)/2;
D = 0.5*p*S*(Cd_i + Cd_min)*(v_lof^2)/2;
T_eff = max_power*e_eff*p_eff_to*sqrt(2)/v_lof;
a = (1/AUW)*(T_eff - D - mu*(AUW*9.81 - L));

to_run = v_lof^2/(2*a)*corr;