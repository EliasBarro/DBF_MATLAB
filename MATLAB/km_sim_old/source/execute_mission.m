function [t_hist, x_hist] = execute_mission(init_states, t_final, rel_tol_perc, abs_tol_vec)
%EXECUTE_MISSION (km_sim toolbox) runs the pre-loaded mission
%
% SYNTAX:
%   [t_hist, x_hist] = execute_mission(init_states, t_final)
%   [t_hist, x_hist] = execute_mission(init_states, t_final, rel_tol_perc, abs_tol_vec)
% 
% INPUTS:
%   init_states   - 14 element vector [x,v,h,aoa,gamma,psi,e,thr,~,~,~,~,north,east]
%                   of initial states. The 4 initial conditions that are
%                   ignored correspond to the termination condition states
%                   which are overriden by the mission plan
%   t_final       - final time in seconds
%   rel_tol_perc  - (optional) gets passed on as the ode45 relative tolerance
%                   default is 1e-5
%   abs_tol_vec   - (optional) gets passed on as the ode45 absolute tolerance 
%                   vector for each state. Defaults is:
%                   [ 1e-1; %x      (m)
%                     1e-2; %v      (m)
%                     1e-2; %h      (m)
%                     2e-5; %aoa   (rad)
%                     1e-3; %gamma (rad)
%                     1e-2; %psi   (rad)
%                     1e-1; %e  (user defined units)
%                     1e-3; %thr    (0-1)
%                     1e-2; %delta st 1
%                     1e-2; %delta st 2
%                     1e-1; %trigger st 1
%                     1e-1; %trigger st 2
%                     1e-1; %north   (m)
%                     1e-1];%east    (m)
% 
% OUTPUTS:
%   t_hist        - time vector that defines when states are output sampled
%   x_hist        - matrix of output states (nt by 14) where nt is the
%                   number of time samples (length of t_hist), columns are 
%                   states.
% 
% see also run_km_sim, mission_phase, current_plane, set_fault

% need to clear end_time so it knows to get new specified end_time
clear phase_to_trans;
clear phase_hold_spd;

% set mission end flag
mission_terminate = false;

% define time history
t_hist      = [];
x_hist      = [];

% sim setup
x_last      = init_states;
t_last      = 0;
run_long    = 10*24*3600;                 % no phase should exceed 10 days
init_dt     = 0.01;
if nargin<3
    rel_tol_perc= 1e-5;
end
if nargin<4
    abs_tol_vec = [ 1e-1; %x
                    1e-2; %v
                    1e-2; %h
                    2e-5; %aoa
                    1e-3; %gamma
                    1e-2; %psi
                    1e-1; %e
                    1e-3; %thr
                    1e-2; %delta st 1
                    1e-2; %delta st 2
                    1e-1; %trigger st 1
                    1e-1; %trigger st 2
                    1e-1; %north
                    1e-1];%east
end
    
sim_opts    = odeset('Events',@transitions,'InitialStep',init_dt,...
    'RelTol',rel_tol_perc,'AbsTol',abs_tol_vec);
[mp_step, mission_plan_out] = mission_phase(); % get current mission_phase
curr_mp     = mission_plan_out{mp_step,1};

% main mission loop
t_st_loop=tic;
while ~mission_terminate 
    % setup terminal transition initializations
    x_last(9)   = -mission_plan_out{mp_step,3}{2};   % delta state 1
    x_last(10)  = -mission_plan_out{mp_step,4}{2};   % delta state 2
    x_last(11) = 1;                                  % reset trigger condition 1
    x_last(12) = 1;                                  % reset trigger condition 2
    
    % run mission_phase
    interval = [t_last, t_last + run_long];
    mp_handle = eval(['@phase_',curr_mp.get_str]);
    [t,x,~,~,indx_event] = ode45(mp_handle, interval, x_last, sim_opts);
    if set_fault('reset'), mission_terminate = true; end     
    
    % if exceeds total mission timer (check #5) then terminate
    if indx_event(1)==5
        mission_terminate = true;
    end
    if mp_step == 23
        wat = 42;
    end
    
    % build up time history
    mp_vec  = double(curr_mp)*ones(length(t),1);
    mps_vec = mp_step*ones(length(t),1);
    x_full  = [x,mp_vec(:),mps_vec(:)];
    t_hist  = [t_hist; t];              %#ok<AGROW>
    x_hist  = [x_hist; x_full];         %#ok<AGROW>

    % setup for next phase
    if ~mission_terminate
        t_last = t(end,:);
        x_last = x(end,:);
        [mp_step, mission_plan_out] = mission_phase(0);
        
        plane = current_plane();
        if any(plane.aero(:,3)<0)
            error(['CD less than zero detected. Please use check_table_data'...
                   'to verify your input script behaves as expected.'])
        end
        
        % execute mission plan sim commands
        flight_phase = false;
        while ~flight_phase
            % get current mission phase (to see if it needs sim processing)
            curr_mp=mission_plan_out{mp_step,1};
            
            % check for goto statement
            if curr_mp == mp.goto
                jump_step = 0;
                try
                    % check delta or trigger state for specific transition
                    if indx_event(1)==1
                        jump_step = eval(mission_plan_out{mp_step,3}); % d1
                    elseif indx_event(1)==2
                        jump_step = eval(mission_plan_out{mp_step,4}); % d2
                    elseif indx_event(1)==3
                        jump_step = eval(mission_plan_out{mp_step,5}); % t1
                    elseif indx_event(1)==4
                        jump_step = eval(mission_plan_out{mp_step,6}); % t2
                    end

                    % default to param state
                    if jump_step==0
                        jump_step = eval(mission_plan_out{mp_step,2});
                    end
                catch
                    error('improperly defined goto mission cmd');
                end

                % process actual jump
                if jump_step == mp_step
                    error('cannot goto current step');
                elseif jump_step == 0
                    error('no goto defined');
                else
                    mp_step = jump_step;
                    [mp_step, mission_plan_out] = mission_phase(mp_step);
                end

            % check for config change
            elseif curr_mp == mp.config
                plane = current_plane(); %#ok<NASGU>
                try
                    eval(mission_plan_out{mp_step,2});
                    eval(mission_plan_out{mp_step,3});
                    eval(mission_plan_out{mp_step,4});
                    eval(mission_plan_out{mp_step,5});
                    eval(mission_plan_out{mp_step,6});
                catch
                    error('improperly defined config change');
                end
                % increment again after the change
                [mp_step, mission_plan_out] = mission_phase(0);
            else
                flight_phase = true;
            end
        end
        
        % check for end of mission plan
        if curr_mp==mp.terminate
            mission_terminate = true; %ya sure I could just "break" but this is cleaner
        end
    end
end
t_elapsed_loop = toc(t_st_loop);
disp(['sim complete, elapsed time: ', num2str(t_elapsed_loop), 'sec']);

% write outputs to base workspace
assignin('base','t_hist',t_hist);
assignin('base','x_hist',x_hist);


% keep t_final scoped for the transition function
function [value,isterminal,direction] = transitions(t,x)  
    value      = [x(9);
                  x(10);
                  x(11);
                  x(12);
                  t_final-t;];
    isterminal = [1;1;1;1;1];
    direction  = [0;0;0;0;0];
end
crit_info(x_hist);
end


