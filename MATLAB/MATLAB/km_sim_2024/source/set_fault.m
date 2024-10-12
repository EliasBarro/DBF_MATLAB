function [fault_out] = set_fault(msg)
%SET_FAULT will display the input "msg" as a warning to the command window
% with the key distinction that this fault will latch. That is, only one
% fault can be set during a sim run, this prevents a warning from being
% output at every time sample of ode45. An unfortunate side effect is that
% additional faults may be shadowed by the first fault. 
% 
% SYNTAX:
%   set fault   -               set_fault(msg)
%   read fault  - [fault_out] = set_fault('')
%   reset fault -               set_fault('reset')
% 
% see also execute_mission, current_plane, mission_phase

persistent fault;

%reset state
if strcmp(msg,'reset')
    fault_out = fault;
    fault = false;  
else
    %if the fault hasn't been set yet, then set the fault
    if ~fault
        fault = true; 
        warning(msg);
    end
    fault_out = fault;
end
    
end