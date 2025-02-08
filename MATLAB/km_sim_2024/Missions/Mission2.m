function Mission5
% THIS IS THE GOOD ONE (12/7/24)
plane = current_plane();
t_final         = 5*60;
init_alt = 668;
turnBackRadius = 250/3.281;
turnAroundRadius = 100/3.281;
laps=1; % CHANGE BASED ON MISSION - IF YOU WANT MAX. # OF LAPS, THEN
         % FIND HOW MANY LAPS CAN BE DONE BEFORE BATTERY STARTS RUNNING
         % OUT. THIS REQUIRES SOME TRIAL & ERROR

% per-plane optional settings (and defaults)
if isfield(plane,'cruise_thr')
    cruise_thr = plane.cruise_thr;
else
    cruise_thr = 1;
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
    normal_load_factor_turn = 2;
end

% Accelerate on the ground at full thrust TODO takeoff speed
mission_plan(1,:) = {mp.ground_thr,  1, {'0',1}, {'0',1}, 'v>6', 'false'}; 

% rotate until weight off wheels (at 2 degrees per second)
mission_plan(2,:) = {mp.ground_rot,  10*pi/180, {'0',1}, {'0',1}, 'g_load>1', 'false'};

% climb at max thr, 20 degree flight path angle, until climbed 10m AGL (delta-h)
mission_plan(3,:) = {mp.climb_max,  20, {'h_dot',10}, {'0',1}, 'false', 'false'};




for i=4: 6:laps*6

% go straight until traveled 500ft from start
mission_plan(i,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east>500/3.28', 'false'}; % 500 feet to meters

% 180 degree right turn
mission_plan(i+1,:) = {mp.turn_radi, -turnBackRadius, {'0',1}, {'psi_dot',-pi}, 'false', 'false'};

% go straight until even with start
mission_plan(i+2,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east<0', 'false'}; 

% 360 degree left turn
mission_plan(i+3,:) = {mp.turn_radi, turnAroundRadius, {'0',1}, {'psi_dot',2 * pi}, 'false', 'false'};

% go straight until 500 ft before start
mission_plan(i+4,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east<-500/3.281', 'false'}; 

% 180 degree right turn
mission_plan(i+5,:) = {mp.turn_radi, -turnBackRadius, {'0',1}, {'psi_dot',-pi}, 'false', 'false'};  % End of final loop


end
%% LANDING

% go straight until even with start
mission_plan(i+6,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east>0', 'false'}; 

% glide at -10 degrees until clearance altitude is less than 1.5m
mission_plan(i+7,:) = {mp.glide, -10, {'0',1}, {'0',1}, ['h<', num2str(init_alt+2)], 'false'};

% glide at -2 degrees until clearance altitude is less than 0.1m
mission_plan(i+8,:) = {mp.glide, -2, {'0',1}, {'0',1}, ['h<', num2str(init_alt+0.2)], 'false'};

% perfect flare right into runway, try to lock into 0 gamma. If getting at max aoa, also quit
mission_plan(i+9,:) = {mp.glide,  0, {'0',1}, {'0',1}, 'gamma>=-1e-4', 'aoa>0.95*max(plane.alt_aero(:,1));'};

% rotate down to level at 2 deg/sec THIS IS CAUSING PLANE TO GO UNDERGROUND
mission_plan(i+10,:) = {mp.ground_rot,  -2*pi/180, {'0',1}, {'0',1}, 'aoa<=0', 'false'};

% roll out
mission_plan(i+11,:) = {mp.ground_thr,  0, {'0',1}, {'0',1}, 'v<5', 'false'};

% terminate
mission_plan(i+12,:) = {mp.terminate,0,0,0,0,0};

% write to mission_phase function
mission_phase(0,mission_plan);

init_states     = zeros(1,14);
init_states(2)  = 0;       % launched at 20m/s
init_states(3)  = init_alt;
init_states(5)  = 0*pi/180; % 10 deg launcher
init_states(6)  = pi/2; % due east
init_states(7)  = plane.e_cap;
init_states(8)  = 1;

% execute
execute_mission(init_states, t_final);