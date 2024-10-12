function [mp_step_out, mission_plan_out] = mission_phase(set_step, mission_in)
%MISSION_PHASE (km_sim toolbox) is mechanism for a pseudo-global variable in 
% the form of a function. Internally, it stores a persistent, settable, 
% current mission plan and step. The mission plan step is "where you are" 
% in the mission plan, specifically, which row you are at.
% 
%	SYNTAX: retrieve the mission_plan  - [mp_step_out, mission_plan_out] = mission_phase();
%           increment the mission phz  -             mission_phase(0);
%           set the mission step       -             mission_phase(desired_step);
%           set the mission            -             mission_phase(0,mission_in);
% 
% see also execute_mission, current_plane, set_fault

% current step within the mission plan
persistent mp_step;

% mission_plan is a cell array where 
% {:,1} are mission phase enumerations
% {:,2} are mission phase parameters or default goto step num
% {:,3:4} are delta state transition conditions
% {:,5:6} are trigger transition conditions
persistent mission_plan; 


%% handle set, reset, increment
if nargin==0
    % retriving step - do nothing (so that else doesn't catch nargin 0)
elseif nargin==1
    if set_step==0 % there is no zero index so this is a flag for just increment
        mp_step = mp_step + 1;
    else
        mp_step = set_step;
    end
elseif nargin==2
    % setting the matrix in the first place
    mission_plan = mission_in;
    mp_step = 1; %init to the first step
else
    error('mission_phase called with incorrect inputs')
end


%% OUTPUTS
mission_plan_out        = mission_plan;
mp_step_out             = mp_step;


end
