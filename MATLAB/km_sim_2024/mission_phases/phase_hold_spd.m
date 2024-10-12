function [states_dot] = phase_hold_spd(t, states)
%phase_hold_spd is called by ode45 to provide the   x_dot = f(x)
% This calculates the power necessary to maintain this speed

% previous proplsn indices
persistent ix;
persistent iy;
persistent iz;
persistent end_time;

if isempty(ix),ix=1;end
if isempty(iy),iy=1;end
if isempty(iz),iz=1;end


%% Input Processing
% retrieve mission parameters
[mp_step, mission_plan_out] = mission_phase();
if isempty(end_time)
    end_time = t + mission_plan_out{mp_step,2};
end

% states are [x, v, h, aoa, gamma, psi, e, thr, delta_tr, trigger_tr];
x       = states(1);    % meters
v       = states(2);	% meters/sec
h       = states(3);    % meters of altitude 
aoa     = states(4);    % radians angle of attack
gamma   = states(5);	% radians of flight path angle
psi     = states(6);	% radians of heading
e       = states(7);	% charge state
thr     = max(min(states(8),1),0);    % 0 - 1 fractional throttle state
north   = states(13);   % north position relative to start
east    = states(14);   %  east position relative to start


%% Define phase dynamics
% get the current plane and AUW
plane   = current_plane();
weight  = plane.empty_W+plane.spec_fuel_W*e;
mass    = weight/9.81;       % N becomes kg

%general parameters
density = model_atmo(h);
q = .5 * density * v^2; % (1/2)*rho*v^2 in Pa
S = plane.S_ref;

% proplsn calcs
if isfield(plane,'thr_lim')
    max_thr = lin_interp1(plane.thr_lim(:,1),plane.thr_lim(:,2),e);
    thr     = min(max_thr, thr);
end
prop                = plane.prop_perf;
try
    [thrust,ix,iy,iz]   = lin_interp3(prop{1},prop{2},prop{3},prop{4},...
        thr,max(0,v),density,ix,iy,iz);
    [e_dot,ix,iy,iz]    = lin_interp3(prop{1},prop{2},prop{3},prop{5},...
        thr,max(0,v),density,ix,iy,iz);
catch
    if v>max(prop{2})
        set_fault(['Velocity out of bounds: vehicle at ' num2str(v) 'm/s but max velocity in propulsion table is ' num2str(max(prop{2})) '. Consider expanding input propulsion table to higher velocities']);
    elseif density<min(prop{3})
        set_fault(['Density out of bounds: vehicle at ' num2str(density) 'kg/m^3, but lowest density in propulsion table is ' num2str(min(prop{3})) '. Consider expanding input propulsion table to lower densities']);
    elseif density>max(prop{3})
        set_fault(['Density out of bounds: vehicle at ' num2str(density) 'kg/m^3, but highest density in propulsion table is ' num2str(max(prop{3})) '. Consider expanding input propulsion table to higher densities']);
    else
        set_fault(['Unpredicted interpolation issue. THR: ' num2str(thr) ',  VEL: ' num2str(v) ', DENS: ' num2str(density)]);
    end
    thrust = 0;
    e_dot = 0;
end

% lift needed maintain current glideslope
lift        = weight*cos(gamma); %ignore thrust sin(aoa) for this: better be smol aoa
g_load      = lift/weight;

% aero calcs (get drag)
CL          = lift/(q*S);
CL = max(min(max(plane.aero(:,2)),CL),min(plane.aero(:,2)));
[CD,inda]   = lin_interp1(plane.aero(:,2),plane.aero(:,3),CL);
drag        = CD*q*S + weight*sin(gamma); % "gravity drag" included
thrust_req  = drag/cos(aoa);
accel       = 0;

% aoa state tracking lags some, only impacts proplsn angle
aoa_target  = lin_interp1(plane.aero(:,2),plane.aero(:,1),CL,inda);
aoa_track_gain = 10;
err_aoa     = aoa-aoa_target;
aoa_dot     = -aoa_track_gain*err_aoa;

% throttle transient (throttle dynamics will lag slightly)
thrust_track_gain  = 0.03*max(max(max(prop{4})));
err_thrust         = thrust-thrust_req;
thr_dot         = -thrust_track_gain*err_thrust;

% Define remaining states_dot
h_dot         = v * sin(gamma);
groundspeed   = v * cos(gamma);
north_dot     = groundspeed*cos(psi);
east_dot      = groundspeed*sin(psi);
psi_dot       = 0;
gamma_dot     = 0;


%% Calculation Monitoring
% end-condition checks
delta_dot1      = eval(mission_plan_out{mp_step,3}{1});
delta_dot2      = eval(mission_plan_out{mp_step,4}{1});
trigger_cond1   = eval(mission_plan_out{mp_step,5});
trigger_cond2   = eval(mission_plan_out{mp_step,6});
if trigger_cond1
    triggr_dot1 = -1e9;
else
    triggr_dot1 = 0;
end
if trigger_cond2
    triggr_dot2 = -1e9;
else
    triggr_dot2 = 0;
end

% standard delta time transition out of this phase
if t>end_time
    triggr_dot1 = -1e9;
end

% sim error checks
if v>75
    error('sim stability issue: diverging speed');
end

% general aircraft faults
if e<0
    set_fault('Vehicle Stored Energy Depleted');
    triggr_dot1 = -1e9;
end


%% assign outputs (states_dot)
states_dot   = [groundspeed;    % 1
                accel;          % 2
                h_dot;          % 3
                aoa_dot;        % 4
                gamma_dot;      % 5
                psi_dot;        % 6
                e_dot;          % 7
                thr_dot;        % 8
                delta_dot1;     % 9
				delta_dot2;     % 10
                triggr_dot1;    % 11
				triggr_dot2;    % 12
                north_dot;      % 13
                east_dot];      % 14
            
if any(isnan(states_dot))
    error('NaN detected in simulation');
end


end