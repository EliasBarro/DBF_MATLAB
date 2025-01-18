function testingMission

plane = current_plane();
t_final         = 20*60;
init_alt = 300;
%turn_radius=30;

% per-plane optional settings (and defaults)
if isfield(plane,'cruise_thr')
    cruise_thr = plane.cruise_thr;
else
    cruise_thr = 0.8;
end
if isfield(plane,'land_thr')
    %land_thr    = plane.land_thr;
else
    %land_thr    = 0.7;
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
    normal_load_factor_turn = 2;
end

% execute takeoff transient for duration specified per-plane. Stop if aoa leveled out or if at 90 percent of total aoa capability
mission_plan(1,:) =  {mp.to_trans, {to_trans_dur, to_trans_gain}, {'0',1}, {'aoa_dot',max(plane.aero(:,1)*to_aoa_frac)},'aoa_dot<0', 'false'};

% climb at max thr, 10 degree flight path angle, until climbed 50m AGL (delta-h)
mission_plan(2,:) = {mp.climb_max,  10, {'h_dot',50}, {'0',1}, 'false', 'false'};

% go straight until traveled 500m from start
mission_plan(3,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'x>500', 'false'};

% 180 degree left turn
mission_plan(4,:) = {mp.turn_gload, normal_load_factor_turn, {'0',1},{'psi_dot',pi},'false','false'};

mission_plan(5,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'x>1200', 'false'};

mission_plan(6,:) = {mp.terminate,0,0,0,0,0};

mission_phase(0,mission_plan);

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