function km()
%KM states are:
%
% x       = states(1);    % meters
% v       = states(2);	  % meters/sec
% h       = states(3);    % meters of altitude 
% aoa     = states(4);    % rad angle of attack
% gamma   = states(5);	  % radians of flight path angle
% psi     = states(6);	  % radians of heading
% e       = states(7);	  % charge state
% thr     = states(8);    % 0 - 1 fractional throttle state
% delta1  = states(9);    % execution control state
% delta2  = states(10);   % execution control state
% triggr1 = states(11);   % execution control state
% triggr2 = states(12);   % execution control state
% north   = states(13);   % north position relative to start
% east    = states(14);   %  east position relative to start
%
% 
% states dot are:
%                 groundspeed    % 1
%                 accel          % 2
%                 h_dot          % 3
%                 aoa_dot        % 4
%                 gamma_dot      % 5
%                 psi_dot        % 6
%                 e_dot          % 7
%                 thr_dot        % 8
%                 delta_dot1     % 9
% 				  delta_dot2     % 10
%                 triggr_dot1    % 11
% 				  triggr_dot2    % 12
%                 north_dot      % 13
%                 east_dot       % 14
%
% these are always defined for every mission phase, otherwise, refer to 
% specific phases for what variables are available for trigger conditions
%
% AVAILABLE PHASES     | mission phase parameter    | thr behavior
%---------------------------------------------------------------------------
%       terminate      |                            | N/A
%       goto           | phase to go to             | N/A
%       config         | statement to run           | N/A
%                      |                            |
%       ground_thr     | thr setting (0-1)          | set
%       ground_rot     | rotate speed (deg/sec)     | hold last thr
%                      |                            |
%       to_trans       | hold time (s), traj_gain   | hold last thr
%       cruise_thr     | thr setting (0-1)          | set
% 		hold_spd       | time to hold (sec)         | dynamically determined
%       climb_max      | gamma (rad)                | max thr
%       glide          | gamma (rad)                |  no thr
%       turn_gload     | g-load of turn             | hold last thr
%       turn_radi      | radius of turn (m)         | hold last thr
%
% see also run_km_sim

end