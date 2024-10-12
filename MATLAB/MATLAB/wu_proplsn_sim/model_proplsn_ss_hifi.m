function [thrust, mech_power, torque_load, angular_rate, current, elec_eff, total_eff] = model_proplsn_ss_hifi(u_fs, thr_const, rho_const, settle_time)
% propulsion model driver, extract last outputs
% note that a corresponding propulsion params file MUST be run prior to
% calling this script to ensure that prop params are loaded in the base
% workspace for use by the propulsion simulink model

%overwrite these params just to be safe that they are correct...
options.LoadExternalInput = 'on';
options.ExternalInput = 'throttle, airspeed_m_s, density_kg_m3';
options.StopTime = 'end_time';

% define input timeseries - (data, time)
assignin('base','airspeed_m_s' ,timeseries([u_fs     , u_fs     ],[0,settle_time]));
assignin('base','throttle'     ,timeseries([thr_const, thr_const],[0,settle_time]));
assignin('base','density_kg_m3',timeseries([rho_const, rho_const],[0,settle_time]));

% define some other basic sim parameters
assignin('base','use_batt' ,true);
assignin('base','use_drain',false);
assignin('base','start_dec',0.8);
assignin('base','sim_dt'   ,0.01);
assignin('base','output_dt',0.25);
assignin('base','end_time' ,settle_time);
modelname = 'standalone_propulsion_model';

% sim
simOut = sim(modelname,options);

% outputs
angular_rate_pre   = simOut.yout{1}.Values;
electric_power_pre = simOut.yout{2}.Values;
thrust_pre         = simOut.yout{3}.Values;
mech_power_pre     = simOut.yout{4}.Values;
current_pre        = simOut.yout{5}.Values;
torque_load_pre    = simOut.yout{6}.Values;

% get only the outputs in question...
temp_ts      = getsampleusingtime(angular_rate_pre, settle_time); 
angular_rate = temp_ts.Data;
temp_ts      = getsampleusingtime(electric_power_pre, settle_time); 
electric_pwr = temp_ts.Data;
temp_ts      = getsampleusingtime(thrust_pre, settle_time); 
thrust       = temp_ts.Data;
temp_ts      = getsampleusingtime(mech_power_pre, settle_time); 
mech_power   = temp_ts.Data;
temp_ts      = getsampleusingtime(current_pre, settle_time); 
current      = temp_ts.Data;
temp_ts      = getsampleusingtime(torque_load_pre, settle_time); 
torque_load  = temp_ts.Data;

% determine total efficiency:
total_eff = (thrust*u_fs)/(electric_pwr);
if total_eff > 1
    total_eff = 1;
elseif total_eff < 0 || isnan(total_eff)
    total_eff = 0;
end

% determine electrical eff:
elec_eff = mech_power/electric_pwr;
if elec_eff > 1
    elec_eff = 1;
elseif elec_eff < 0 || isnan(elec_eff)
    elec_eff = 0;
end

end