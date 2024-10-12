function states_dot = phase_ground_thr(t, states)              %#ok<INUSL>
% phase_takeoff is called by ode45 to provide the   x_dot = f(x)
% This calculates the acceleration along the path for a takeoff roll.

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
thr_command                 = mission_plan_out{mp_step,2};

% states are [x, v, h, aoa, gamma, psi, e, thr, delta_tr, trigger_tr];
x       = states(1);    % meters
v       = states(2);	% meters/sec
h       = states(3);    % meters of altitude 
aoa     = states(4);    % rad angle of attack
gamma   = states(5);	% radians of flight path angle
psi     = states(6);	% radians of heading
e       = states(7);	% charge state
thr     = states(8);    % 0 - 1 fractional throttle state
north   = states(13);   % north position relative to start
east    = states(14);   %  east position relative to start


%% Define Dynamics for this phase
% throttle transient
thr_track_gain  = 10;
err_thr         = thr-thr_command;
thr_dot         = -thr_track_gain*err_thr;
thr             = min(max(0,thr),1); % limit the used throttle to 0-1

% get the current plane and AUW
plane   = current_plane();
weight  = plane.empty_W+plane.spec_fuel_W*e;
mass    = weight/9.81;       % N becomes kg

%general parameters
density = model_atmo(h);
q       = 0.5 * density * v^2;    % (1/2)*rho*v^2 in Pa
S       = plane.S_ref;            % so obviously S_ref needs to be in meters squared
max_aoa = max(plane.aero(:,1));
min_aoa = min(plane.aero(:,1));

% proplsn calcs
if isfield(plane,'thr_lim')
    max_thr = lin_interp1(plane.thr_lim(:,1),plane.thr_lim(:,2),e);
    thr     = min(max_thr, thr);
end
prop                = plane.prop_perf;
[thrust,ix,iy,iz]   = lin_interp3(prop{1},prop{2},prop{3},prop{4},...
    thr,max(0,v),density,ix,iy,iz);
[e_dot, ix,iy,iz]   = lin_interp3(prop{1},prop{2},prop{3},prop{5},...
    thr,max(0,v),density,ix,iy,iz);

% aero calcs
aero_aoa    = min(max(min_aoa,aoa),max_aoa); % aoa limiting
[CL,inda]   = lin_interp1(plane.aero(:,1),plane.aero(:,2),aero_aoa);
CD          = lin_interp1(plane.aero(:,1),plane.aero(:,3),aero_aoa,inda);
normal_force= CL*q*S + thrust*sin(aoa)-weight*cos(gamma); %total including thrust
drag        = CD*q*S;
g_load      = CL*q*S/weight;

% friction calc
friction = plane.TO_frict_coeff*abs(normal_force);
if abs(v)<1e-9 && thrust*cos(aoa)<friction
    friction = thrust*cos(aoa)-weight*sin(gamma);
end

% Define remaining states dot
accel           = (thrust*cos(aoa) - weight*sin(gamma) - drag - friction)/mass; % m/s^2
groundspeed     = v * cos(gamma);
h_dot           = v * sin(gamma);
gamma_dot       = 0;
aoa_dot         = 0;
psi_dot         = 0;
north_dot       = groundspeed*cos(psi);
east_dot        = groundspeed*sin(psi);


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

% phase specific faults:
if normal_force>0 % check for takeoff
    triggr_dot1 = -1e9; % advance to something hopefully...
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
if isnan(states_dot)
    error('NaN detected in simulation');
end


end