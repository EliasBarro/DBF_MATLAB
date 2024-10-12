function Mission1

plane = current_plane();
t_final         = 10*60;
init_alt = 669;
turnBackRadius = 250/3.281;
turnAroundRadius = 100/3.281;

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

% go straight until traveled 500ft from start
mission_plan(4,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'x>500/3.28', 'false'}; % 500 feet to meters

% 180 degree right turn
mission_plan(5,:) = {mp.turn_radi, -turnBackRadius, {'0',1}, {'psi_dot',-pi}, 'false', 'false'};

% go straight until even with start
mission_plan(6,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east<0', 'false'}; 

% 360 degree left turn
mission_plan(7,:) = {mp.turn_radi, turnAroundRadius, {'0',1}, {'psi_dot',2 * pi}, 'false', 'false'};

% go straight until even with start
mission_plan(8,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east<-500/3.281', 'false'}; 

% 180 degree right turn
mission_plan(9,:) = {mp.turn_radi, -turnBackRadius, {'0',1}, {'psi_dot',-pi}, 'false', 'false'};

% go straight until 500 from start
mission_plan(10,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east>500/3.281', 'false'}; % End of first loop


% 180 degree right turn
mission_plan(11,:) = {mp.turn_radi, -turnBackRadius, {'0',1}, {'psi_dot',-pi}, 'false', 'false'};

% go straight until even with start
mission_plan(12,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east<0', 'false'}; 

% 360 degree left turn
mission_plan(13,:) = {mp.turn_radi, turnAroundRadius, {'0',1}, {'psi_dot',2 * pi}, 'false', 'false'};

% go straight until even with start
mission_plan(14,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east<-500/3.281', 'false'}; 

% 180 degree right turn
mission_plan(15,:) = {mp.turn_radi, -turnBackRadius, {'0',1}, {'psi_dot',-pi}, 'false', 'false'};

% go straight until 500 from start
mission_plan(16,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east>500/3.281', 'false'}; % End of second loop


% 180 degree right turn
mission_plan(17,:) = {mp.turn_radi, -turnBackRadius, {'0',1}, {'psi_dot',-pi}, 'false', 'false'};

% go straight until even with start
mission_plan(18,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east<0', 'false'}; 

% 360 degree left turn
mission_plan(19,:) = {mp.turn_radi, turnAroundRadius, {'0',1}, {'psi_dot',2 * pi}, 'false', 'false'};

% go straight until even with start
mission_plan(20,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east<-500/3.281', 'false'}; 

% 180 degree right turn
mission_plan(21,:) = {mp.turn_radi, -turnBackRadius, {'0',1}, {'psi_dot',-pi}, 'false', 'false'};  % End of final loop

% go straight until 500 from start
mission_plan(22,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east>500/3.281', 'false'}; % End of second loop


% 180 degree right turn
mission_plan(23,:) = {mp.turn_radi, -turnBackRadius, {'0',1}, {'psi_dot',-pi}, 'false', 'false'};

% go straight until even with start
mission_plan(24,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east<0', 'false'}; 

% 360 degree left turn
mission_plan(25,:) = {mp.turn_radi, turnAroundRadius, {'0',1}, {'psi_dot',2 * pi}, 'false', 'false'};

% go straight until even with start
mission_plan(26,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east<-500/3.281', 'false'}; 

% 180 degree right turn
mission_plan(27,:) = {mp.turn_radi, -turnBackRadius, {'0',1}, {'psi_dot',-pi}, 'false', 'false'};  % End of final loop

% go straight until 500 from start
mission_plan(28,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east>500/3.281', 'false'}; % End of second loop


% 180 degree right turn
mission_plan(29,:) = {mp.turn_radi, -turnBackRadius, {'0',1}, {'psi_dot',-pi}, 'false', 'false'};

% go straight until even with start
mission_plan(30,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east<0', 'false'}; 

% 360 degree left turn
mission_plan(31,:) = {mp.turn_radi, turnAroundRadius, {'0',1}, {'psi_dot',2 * pi}, 'false', 'false'};

% go straight until even with start
mission_plan(32,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east<-500/3.281', 'false'}; 

% 180 degree right turn
mission_plan(33,:) = {mp.turn_radi, -turnBackRadius, {'0',1}, {'psi_dot',-pi}, 'false', 'false'};  % End of final loop

% go straight until 500 from start
mission_plan(34,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east>500/3.281', 'false'}; % End of second loop


% 180 degree right turn
mission_plan(35,:) = {mp.turn_radi, -turnBackRadius, {'0',1}, {'psi_dot',-pi}, 'false', 'false'};

% go straight until even with start
mission_plan(36,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east<0', 'false'}; 

% 360 degree left turn
mission_plan(37,:) = {mp.turn_radi, turnAroundRadius, {'0',1}, {'psi_dot',2 * pi}, 'false', 'false'};

% go straight until even with start
mission_plan(38,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east<-500/3.281', 'false'}; 

% 180 degree right turn
mission_plan(39,:) = {mp.turn_radi, -turnBackRadius, {'0',1}, {'psi_dot',-pi}, 'false', 'false'};  % End of final loop

% go straight until 500 from start
mission_plan(40,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east>500/3.281', 'false'}; % End of second loop


% 180 degree right turn
mission_plan(41,:) = {mp.turn_radi, -turnBackRadius, {'0',1}, {'psi_dot',-pi}, 'false', 'false'};

% go straight until even with start
mission_plan(42,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east<0', 'false'}; 

% 360 degree left turn
mission_plan(43,:) = {mp.turn_radi, turnAroundRadius, {'0',1}, {'psi_dot',2 * pi}, 'false', 'false'};

% go straight until even with start
mission_plan(44,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east<-500/3.281', 'false'}; 

% 180 degree right turn
mission_plan(45,:) = {mp.turn_radi, -turnBackRadius, {'0',1}, {'psi_dot',-pi}, 'false', 'false'};  % End of final loop

% go straight until 500 from start
mission_plan(46,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east>500/3.281', 'false'}; % End of second loop


% 180 degree right turn
mission_plan(47,:) = {mp.turn_radi, -turnBackRadius, {'0',1}, {'psi_dot',-pi}, 'false', 'false'};

% go straight until even with start
mission_plan(48,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east<0', 'false'}; 

% 360 degree left turn
mission_plan(49,:) = {mp.turn_radi, turnAroundRadius, {'0',1}, {'psi_dot',2 * pi}, 'false', 'false'};

% go straight until even with start
mission_plan(50,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east<-500/3.281', 'false'}; 

% 180 degree right turn
mission_plan(51,:) = {mp.turn_radi, -turnBackRadius, {'0',1}, {'psi_dot',-pi}, 'false', 'false'};  % End of final loop

% go straight until 500 from start
mission_plan(52,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east>500/3.281', 'false'}; % End of second loop


% 180 degree right turn
mission_plan(53,:) = {mp.turn_radi, -turnBackRadius, {'0',1}, {'psi_dot',-pi}, 'false', 'false'};

% go straight until even with start
mission_plan(54,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east<0', 'false'}; 

% 360 degree left turn
mission_plan(55,:) = {mp.turn_radi, turnAroundRadius, {'0',1}, {'psi_dot',2 * pi}, 'false', 'false'};

% go straight until even with start
mission_plan(56,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east<-500/3.281', 'false'}; 

% 180 degree right turn
mission_plan(57,:) = {mp.turn_radi, -turnBackRadius, {'0',1}, {'psi_dot',-pi}, 'false', 'false'};  % End of final loop

% go straight until 500 from start
mission_plan(58,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east>500/3.281', 'false'}; % End of second loop


% 180 degree right turn
mission_plan(59,:) = {mp.turn_radi, -turnBackRadius, {'0',1}, {'psi_dot',-pi}, 'false', 'false'};

% go straight until even with start
mission_plan(60,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east<0', 'false'}; 

% 360 degree left turn
mission_plan(61,:) = {mp.turn_radi, turnAroundRadius, {'0',1}, {'psi_dot',2 * pi}, 'false', 'false'};

% go straight until even with start
mission_plan(62,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east<-500/3.281', 'false'}; 

% 180 degree right turn
mission_plan(63,:) = {mp.turn_radi, -turnBackRadius, {'0',1}, {'psi_dot',-pi}, 'false', 'false'};  % End of final loop

% go straight until 500 from start
mission_plan(64,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east>500/3.281', 'false'}; % End of second loop


% 180 degree right turn
mission_plan(65,:) = {mp.turn_radi, -turnBackRadius, {'0',1}, {'psi_dot',-pi}, 'false', 'false'};

% go straight until even with start
mission_plan(66,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east<0', 'false'}; 

% 360 degree left turn
mission_plan(67,:) = {mp.turn_radi, turnAroundRadius, {'0',1}, {'psi_dot',2 * pi}, 'false', 'false'};

% go straight until even with start
mission_plan(68,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east<-500/3.281', 'false'}; 

% 180 degree right turn
mission_plan(69,:) = {mp.turn_radi, -turnBackRadius, {'0',1}, {'psi_dot',-pi}, 'false', 'false'};  % End of final loop

% go straight until 500 from start
mission_plan(70,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east>500/3.281', 'false'}; % End of second loop


% 180 degree right turn
mission_plan(71,:) = {mp.turn_radi, -turnBackRadius, {'0',1}, {'psi_dot',-pi}, 'false', 'false'};

% go straight until even with start
mission_plan(72,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east<0', 'false'}; 

% 360 degree left turn
mission_plan(73,:) = {mp.turn_radi, turnAroundRadius, {'0',1}, {'psi_dot',2 * pi}, 'false', 'false'};

% go straight until even with start
mission_plan(74,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east<-500/3.281', 'false'}; 

% 180 degree right turn
mission_plan(75,:) = {mp.turn_radi, -turnBackRadius, {'0',1}, {'psi_dot',-pi}, 'false', 'false'};  % End of final loop

% go straight until 500 from start
mission_plan(76,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east>500/3.281', 'false'}; % End of second loop


% 180 degree right turn
mission_plan(77,:) = {mp.turn_radi, -turnBackRadius, {'0',1}, {'psi_dot',-pi}, 'false', 'false'};

% go straight until even with start
mission_plan(78,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east<0', 'false'}; 

% 360 degree left turn
mission_plan(79,:) = {mp.turn_radi, turnAroundRadius, {'0',1}, {'psi_dot',2 * pi}, 'false', 'false'};

% go straight until even with start
mission_plan(80,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east<-500/3.281', 'false'}; 

% 180 degree right turn
mission_plan(81,:) = {mp.turn_radi, -turnBackRadius, {'0',1}, {'psi_dot',-pi}, 'false', 'false'};  % End of final loop

% go straight until 500 from start
mission_plan(82,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east>500/3.281', 'false'}; % End of second loop


% 180 degree right turn
mission_plan(83,:) = {mp.turn_radi, -turnBackRadius, {'0',1}, {'psi_dot',-pi}, 'false', 'false'};

% go straight until even with start
mission_plan(84,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east<0', 'false'}; 

% 360 degree left turn
mission_plan(85,:) = {mp.turn_radi, turnAroundRadius, {'0',1}, {'psi_dot',2 * pi}, 'false', 'false'};

% go straight until even with start
mission_plan(86,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east<-500/3.281', 'false'}; 

% 180 degree right turn
mission_plan(87,:) = {mp.turn_radi, -turnBackRadius, {'0',1}, {'psi_dot',-pi}, 'false', 'false'};  % End of final loop

% go straight until 500 from start
mission_plan(88,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east>500/3.281', 'false'}; % End of second loop


% 180 degree right turn
mission_plan(89,:) = {mp.turn_radi, -turnBackRadius, {'0',1}, {'psi_dot',-pi}, 'false', 'false'};

% go straight until even with start
mission_plan(90,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east<0', 'false'}; 

% 360 degree left turn
mission_plan(91,:) = {mp.turn_radi, turnAroundRadius, {'0',1}, {'psi_dot',2 * pi}, 'false', 'false'};

% go straight until even with start
mission_plan(92,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east<-500/3.281', 'false'}; 

% 180 degree right turn
mission_plan(93,:) = {mp.turn_radi, -turnBackRadius, {'0',1}, {'psi_dot',-pi}, 'false', 'false'};  % End of final loop

% go straight until 500 from start
mission_plan(94,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east>500/3.281', 'false'}; % End of second loop


% 180 degree right turn
mission_plan(95,:) = {mp.turn_radi, -turnBackRadius, {'0',1}, {'psi_dot',-pi}, 'false', 'false'};

% go straight until even with start
mission_plan(96,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east<0', 'false'}; 

% 360 degree left turn
mission_plan(97,:) = {mp.turn_radi, turnAroundRadius, {'0',1}, {'psi_dot',2 * pi}, 'false', 'false'};

% go straight until even with start
mission_plan(98,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east<-500/3.281', 'false'}; 

% 180 degree right turn
mission_plan(99,:) = {mp.turn_radi, -turnBackRadius, {'0',1}, {'psi_dot',-pi}, 'false', 'false'};  % End of final loop

% go straight until 500 from start
mission_plan(100,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east>500/3.281', 'false'}; % End of second loop


% 180 degree right turn
mission_plan(101,:) = {mp.turn_radi, -turnBackRadius, {'0',1}, {'psi_dot',-pi}, 'false', 'false'};

% go straight until even with start
mission_plan(102,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east<0', 'false'}; 

% 360 degree left turn
mission_plan(103,:) = {mp.turn_radi, turnAroundRadius, {'0',1}, {'psi_dot',2 * pi}, 'false', 'false'};

% go straight until even with start
mission_plan(104,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east<-500/3.281', 'false'}; 

% 180 degree right turn
mission_plan(105,:) = {mp.turn_radi, -turnBackRadius, {'0',1}, {'psi_dot',-pi}, 'false', 'false'};  % End of final loop

% go straight until even with start
mission_plan(106,:) = {mp.cruise_thr, cruise_thr, {'0',1}, {'0',1}, 'east<0', 'false'}; 
%% LANDING
% glide at -10 degrees until clearance altitude is less than 1.5m
mission_plan(107,:) = {mp.glide, -10, {'0',1}, {'0',1}, ['h<', num2str(init_alt+2)], 'false'};

% glide at -2 degrees until clearance altitude is less than 0.1m
mission_plan(108,:) = {mp.glide, -2, {'0',1}, {'0',1}, ['h<', num2str(init_alt+0.2)], 'false'};

% perfect flare right into runway, try to lock into 0 gamma. If getting at max aoa, also quit
mission_plan(109,:) = {mp.glide,  0, {'0',1}, {'0',1}, 'gamma>=-1e-4', 'aoa>0.95*max(plane.alt_aero(:,1));'};

% rotate down to level at 2 deg/sec
mission_plan(110,:) = {mp.ground_rot,  -2*pi/180, {'0',1}, {'0',1}, 'aoa<=0', 'false'};

% roll out
mission_plan(111,:) = {mp.ground_thr,  0, {'0',1}, {'0',1}, 'v<5', 'false'};

% terminate
mission_plan(112,:) = {mp.terminate,0,0,0,0,0};

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