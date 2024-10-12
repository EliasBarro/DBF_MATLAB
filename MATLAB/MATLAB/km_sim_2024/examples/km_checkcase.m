function km_checkcase
%KM_CHECKCASE runs through all of the km_sim capabilities in a single
%mission
%
% note: turns to the right are positive

% constants
plane           = current_plane();
t_final         = 20*60;
init_alt        = 300;    % start at 300m density altitude
turn_radius     = 30;     % meters

% per-plane optional settings (and defaults)
if isfield(plane,'cruise_thr')
    cruise_thr = plane.cruise_thr;
else
    cruise_thr = 0.8;
end
if isfield(plane,'land_thr')
    land_thr    = plane.land_thr;
else
    land_thr    = 0.7;
end
if isfield(plane,'to_trans_dur')
    to_trans_dur = plane.to_trans_dur;
else
    to_trans_dur = 1;
end
if isfield(plane,'to_trans_gain')
    to_trans_gain = plane.to_trans_gain;
else
    to_trans_gain = 2;
end
if isfield(plane,'to_aoa_frac')
    to_aoa_frac = plane.to_aoa_frac;
else
    to_aoa_frac = 0.9;
end
if isfield(plane,'maneuver_g')
    normal_load_factor_turn = plane.maneuver_g;
else
    normal_load_factor_turn = 2.5;
end


%% LAUNCHER TAKEOFF
% format is: {mp.phase, phase parameter, delta state condition 1{state,delta}, delta state condition 2{state,delta}, trigger condition 1, trigger condition 2}

% execute takeoff transient for duration specified per-plane. Stop if aoa leveled out or if at 90 percent of total aoa capability
mission_plan(1,:) =  {mp.to_trans, {to_trans_dur, to_trans_gain}, {'0',1}, {'aoa_dot',max(plane.aero(:,1)*to_aoa_frac)},'aoa_dot<0', 'false'};

% climb at max thr, 10 degree flight path angle, until climbed 50m AGL (delta-h)
mission_plan(2,:) = {mp.climb_max,  10, {'h_dot',50}, {'0',1}, 'false', 'false'};

% config change - remove flaps (alternate aero is cruise config)
mission_plan(3,:) = {mp.config, 'set_alt_aero;','','','',''};

% go straight until traveled 500m from start
mission_plan(4,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'x>500', 'false'};

% 90 degree left turn
mission_plan(5,:) = {mp.turn_gload, -normal_load_factor_turn, {'0',1},{'psi_dot',-pi/2},'false','false'};

% go straight until at 1000m North
mission_plan(6,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'north>500', 'false'};


%% LOITER in racetrack pattern until 3 minutes of flight
% constant g left turn sweep through 180degrees of heading change
mission_plan(7,:) = {mp.turn_gload, -normal_load_factor_turn, {'0',1}, {'psi_dot',-pi}, 'false', 'false'};

% cruise for 50m
mission_plan(8,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'groundspeed',50}, 't>3*60', 'false'};

% CHECK1 - if exceeded 3 min, then skip to landing approach, else proceed normally
mission_plan(9,:) = {mp.goto, '10','0','0','16','0'};

% cruise for 100m
mission_plan(10,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'groundspeed',100}, 'false', 'false'};

% constant g left turn sweep through 180 degrees of heading change
mission_plan(11,:) = {mp.turn_gload, -normal_load_factor_turn, {'0',1}, {'psi_dot',-pi}, 'false', 'false'};

% cruise for 50m
mission_plan(12,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'groundspeed',50}, 't>3*60', 'false'};

% CHECK 2 - if exceeded 3 mins, then land, else proceed normally
mission_plan(13,:) = {mp.goto,'14','0','0','18','0'};

% cruise until at 500m North
mission_plan(14,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'north>500', 'false'};

% goto indx 7 - the first turn
mission_plan(15,:) = {mp.goto, '7','0','0','0','0'};

% with the above goto we have created no way to naturally transition
% to the next segment, this makes it the natural place to put our CHECK1&2
% actions


%% RUNWAY LINEUP
% southbound (CHECK1) right 90deg turn to head west
mission_plan(16,:) = {mp.turn_gload, normal_load_factor_turn, {'0',1}, {'psi_dot',pi/2}, 'false', 'false'};
mission_plan(17,:) = {mp.goto, '19','0','0','0','0'};

% northbound (CHECK2) left 90deg turn to head west
mission_plan(18,:) = {mp.turn_gload, -normal_load_factor_turn, {'0',1}, {'psi_dot',-pi/2}, 'false', 'false'};

% fly straight until at 250m west of runway
mission_plan(19,:) = {mp.cruise_thr, land_thr, {'0',1}, {'0',1}, 'east<-250', 'false'};

% turn south (left 90deg turn)
mission_plan(20,:) = {mp.turn_gload, -normal_load_factor_turn, {'0',1}, {'psi_dot',-pi/2}, 'false', 'false'};

% fly straight, setup for final landing turn
mission_plan(21,:) = {mp.cruise_thr, land_thr, {'0',1}, {'0',1}, ['north<',num2str(turn_radius)], 'false'};

% constant-radius left turn to line up with runway
mission_plan(22,:) = {mp.turn_radi, -turn_radius, {'0',1}, {'psi_dot',-pi/2}, 'false', 'false'};

% fly until at descent entry point
mission_plan(23,:) = {mp.cruise_thr, land_thr, {'0',1}, {'0',1}, 'east>-150', 'false'};


%% LANDING
% glide at -10 degrees until clearance altitude is less than 1.5m
mission_plan(24,:) = {mp.glide, -10, {'0',1}, {'0',1}, ['h<', num2str(init_alt+2)], 'false'};

% glide at -2 degrees until clearance altitude is less than 0.1m
mission_plan(25,:) = {mp.glide, -2, {'0',1}, {'0',1}, ['h<', num2str(init_alt+0.2)], 'false'};

% perfect flare right into runway, try to lock into 0 gamma. If getting at max aoa, also quit
mission_plan(26,:) = {mp.glide,  0, {'0',1}, {'0',1}, 'gamma>=-1e-4', 'aoa>0.95*max(plane.aero(:,1));'};

% rotate down to level at 2 deg/sec
mission_plan(27,:) = {mp.ground_rot,  -2*pi/180, {'0',1}, {'0',1}, 'aoa<=0', 'false'};

% roll out
mission_plan(28,:) = {mp.ground_thr,  0, {'0',1}, {'0',1}, 'v<0.1', 'false'};


%% SECOND TAKEOFF (ROLLING)
% accelerate to 25m/s
mission_plan(29,:) = {mp.ground_thr,  1, {'0',1}, {'0',1}, 'v>25', 'false'};

% rotate until weight off wheels (at 2 degrees per second)
mission_plan(30,:) = {mp.ground_rot,  2*pi/180, {'0',1}, {'0',1}, 'g_load>1', 'false'};

% climb at max thr, 10 degree flight path angle, until climbed 100m AGL (delta-h)
mission_plan(31,:) = {mp.climb_max,  10, {'h_dot',100}, {'0',1}, 'false', 'false'};

% turn around
mission_plan(32,:) = {mp.turn_gload, -normal_load_factor_turn, {'0',1}, {'psi_dot',-pi}, 'false', 'false'};


%% CRUISE AT 25m/s for 2 minutes
% dummy "cruise" phase just to check if we are above or below 25m/s
% GOTO a zero thr or full throttle accordingly
mission_plan(33,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'groundspeed',50}, 'v<25', 'v>25'};
mission_plan(34,:) = {mp.goto,'38','0','0','35','37'};

% accelerate
mission_plan(35,:) = {mp.cruise_thr, 1, {'0',1}, {'0',1}, 'v>=25', 'false'};
mission_plan(36,:) = {mp.goto,'38','0','0','0','0'};

% decelerate
mission_plan(37,:) = {mp.cruise_thr, 0, {'0',1}, {'0',1}, 'v<=25', 'false'};

% hold speed for 3 minutes
mission_plan(38,:) = {mp.hold_spd, 2*60, {'0',1}, {'0',1}, 'false', 'false'};

% if holding again, use the following line to reset the timer:
% mission_plan(<next>,:) = {mp.config,'clear phase_hold_spd;','','','',''};

% turn around
mission_plan(39,:) = {mp.turn_gload, normal_load_factor_turn, {'0',1}, {'psi_dot',pi}, 'false', 'false'};


%% STOP
% terminate
mission_plan(40,:) = {mp.terminate,0,0,0,0,0};

% write to mission_phase function
mission_phase(0,mission_plan);


%% initialize and execute mission
% states are [x, v, h, aoa, gamma, psi, e, thr, delta_1, delta_2, trigr_1, trigr_2, north, east];
% x       = states(1);  % meters
% v       = states(2);	% meters/sec
% h       = states(3);  % meters of altitude 
% aoa     = states(4);  % radians angle of attack
% gamma   = states(5);	% radians of flight path angle
% psi     = states(6);	% radians of heading
% e       = states(7);	% charge state
% thr     = states(8);  % 0 - 1 fractional throttle state
% delta1  = states(9);
% delta2  = states(10);
% trigr1  = states(11);
% trigr2  = states(12);
% north   = states(13);
% east    = states(14);

init_states     = zeros(1,14);
init_states(2)  = 15;       % launched at 15m/s
init_states(3)  = init_alt;
init_states(5)  = 10*pi/180; % 10 deg launcher
init_states(6)  = pi/2; % due east
init_states(7)  = plane.e_cap;
init_states(8)  = 1;

% execute
execute_mission(init_states, t_final);


end