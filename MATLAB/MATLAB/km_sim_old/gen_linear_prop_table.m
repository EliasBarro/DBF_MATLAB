function prop_perf = gen_linear_prop_table(sl_static_edot, max_alt, thrust_vs_velo)
%GEN_LINEAR_PROP_TABLE generates a simple linear thrust table for km_sim
%(modified)
%
% SYNTAX:
%   gen_linear_prop_table(sl_static_edot, max_alt, thrust_vs_velo)
%
% INPUTS:
%   sl_static_edot      -   The sea-level static energy_dot of the power system
%                           in units consistent with the energy units of
%                           the project (for mAh energy units, take 
%                           Amps*1000/3600 to get mAh per second).
%   max_alt             -   Maximum altitude the table should be valid for
%                           in meters
%   thrust_vs_velo      -   CSV file containing thrust values in the first
%                           column and corresponding velocity values in the 
%                           second column
%
% OUTPUTS:
%   prop_perf   - A cell array containing a 3D thrust table, a 3D "eDot"
%                 table, and three vectors of breakpoints to define the
%                 elements for each table
%       prop_perf{1} = bp_thr   -   breakpoints for dimension 1, throttle
%       prop_perf{2} = bp_u_fs  -   breakpoitns for dim 2, forward speed (m/s)
%       prop_perf{3} = bp_dens  -   breakpoints for dim 3, air density(kg/m^3)
%       prop_perf{4} = thrust_3D
%       prop_perf{5} = eDot_3D
%
% for <1000m, the table will only have 3 density points. For above 1000m, the table
% will have 10 points.


%% SETUP
% generate density breakpoints
if max_alt<1000
    alt_v = [0 500 1000];
else
    alt_v = linspace(0,max_alt,6);  % change
end
bp_dens = zeros(1,length(alt_v));
for indx=1:length(alt_v)
    bp_dens(indx) = model_atmo(alt_v(indx));
end

bp_dens = fliplr(bp_dens); % change

% generate other breakpoints
% import velo_vs_thrust data
th_v_velo = importdata(thrust_vs_velo);
th_v_velo_hold = th_v_velo.data;
veloimport = th_v_velo_hold(:,2);
bp_u_fs = rmmissing(veloimport);

thrustimport = th_v_velo_hold(:,1);
thrust = rmmissing(thrustimport);

bp_thr  = [0 0.2 0.4 0.6 0.8 1];
%% GENERATE TABLES
% thrust table
thrust_3D   = zeros(length(bp_thr), length(bp_u_fs), length(bp_dens));

for indz=1:length(bp_dens)
    this_dens = bp_dens(indz);
    for indy=1:length(bp_u_fs)
        for indx=1:length(bp_thr)
            this_thr = bp_thr(indx);
            thrust_3D(indx,indy,indz) = 9.81*thrust(indy) ...
                *(this_thr/1)*(this_dens/1.225);
        end
    end
end

% eDot (power) table
eDot_3D     = thrust_3D.*sl_static_edot/(9.81*thrust(1));


%% FINAL OUTPUT
prop_perf = {   bp_thr, bp_u_fs, bp_dens,... % breakpoints first
                    thrust_3D, eDot_3D};     % data tables second
                
                
end