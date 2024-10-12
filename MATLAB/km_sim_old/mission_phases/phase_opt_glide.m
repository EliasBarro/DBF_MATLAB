function [states_dot] = phase_opt_glide(t, states)             %#ok<INUSL>
% THIS HASN"T BEEN TOUCHED YET
% This calculates the acceleration along the path for a climb at max thr

% previous proplsn indices
persistent ix;
persistent iy;
persistent iz;

if isempty(ix),ix=1;end
if isempty(iy),iy=1;end
if isempty(iz),iz=1;end


%% Input Processing
% retrieve mission parameters
[mp_step, mission_plan_out] = mission_phase();
v_target                    = mission_plan_out{mp_step,2};

% states are [x, v, h, aoa, gamma, psi, e, thr, delta_tr, trigger_tr];
x       = states(1);    % meters
v       = states(2);	% meters/sec
h       = states(3);    % meters of altitude 
aoa     = states(4);    % radians angle of attack
gamma   = states(5);	% radians of flight path angle
psi     = states(6);	% radians of heading
e       = states(7);	% charge state
thr     = states(8);    % 0 - 1 fractional throttle state
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

% throttle transient
thr_track_gain  = 10;
err_thr         = thr-1; % max throttle (1) is hardcoded into climb
thr_dot         = -thr_track_gain*err_thr;
thr             = min(max(0,thr),1); % limit the used throttle to 0-1

% proplsn calcs
if isfield(plane,'thr_lim')
    max_thr = lin_interp1(plane.thr_lim(:,1),plane.thr_lim(:,2),e);
    thr     = min(max_thr, thr);
end
prop            = plane.prop_perf;
[thrust,ix,iy,iz]   = lin_interp3(prop{1},prop{2},prop{3},prop{4},...
    thr,max(0,v),density,ix,iy,iz);
[e_dot,ix,iy,iz]    = lin_interp3(prop{1},prop{2},prop{3},prop{5},...
    thr,max(0,v),density,ix,iy,iz);

% gamma_dot to determine the extra lift for changing gamma at the start
gamma_track_gain    = 0.1;
err_v               = v-v_target;
gamma_dot           = gamma_track_gain*err_v;
lift                = mass*v*gamma_dot; %F=ma - gamma dot needs to be radps

% add lift needed maintain current glideslope
lift = lift + weight*cos(gamma)-thrust*sin(aoa);

% aero calcs (get drag)
CL          = lift/(q*S);
[CD,inda]   = lin_interp1(plane.aero(:,2),plane.aero(:,3),CL);
drag        = CD*q*S - weight*sin(gamma); % "gravity drag" included
g_load      = lift/weight;

% aoa state tracking lags some, only impacts proplsn angle
aoa_target  = lin_interp1(plane.aero(:,2),plane.aero(:,1),CL,inda);
aoa_track_gain = 10;
err_aoa     = aoa-aoa_target;
aoa_dot     = -aoa_track_gain*err_aoa;

% Define remaining states_dot
accel         = (thrust*cos(aoa) - drag)/mass; % m/s^2
h_dot         = v * sin(gamma);
groundspeed   = v * cos(gamma);
psi_dot       = 0;

%% Calculation Monitoring
% end-condition checks
delta_dot       = eval(mission_plan_out{mp_step,3}{1});
if eval(mission_plan_out{mp_step,4})
    triggr_dot1 = -1e9;
else
    triggr_dot1 = 0;
end

% phase specific faults:
if v<1
    set_fault('vehicle cannot climb');
    triggr_dot1 = -1e9;
end

% sim error checks
if v>75
    error('sim stability issue: diverging speed');
end

% general aircraft faults
max_lift = max(plane.aero(:,2))*q*S;
if (lift>max_lift)
    set_fault('vehicle cannot produce the lift needed');
    triggr_dot1 = -1e9;
elseif e<0
    set_fault('Vehicle Stored Energy Depleted');
    triggr_dot1 = -1e9;
elseif g_load > plane.g_limit
    set_fault('g limit exceeded, vehicle failed structurally');
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
                delta_dot;      % 9
                triggr_dot1;    % 10
                north_dot;      % 11
                east_dot];      % 12
if isnan(states_dot)
    error('NaN detected in simulation');
end


end