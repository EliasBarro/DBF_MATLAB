function run_km_sim(vehicle_file, mission_file, post_pro_file)
%RUN_KM_SIM is the top-level function for executing km_sim, which stands
%for "kinematic mission simulation".
%
% SYNTAX:
%   run_km_sim(vehicle_file, mission_file)
%   run_km_sim(vehicle_file, mission_file, post_pro_file)
%
% INPUTS:
%   vehicle_file -  this is a .mat file that contains the a struct variable
%                   named "plane" which specifies the aircraft tables and data
%   mission_file -  this is a function that defines the mission plan and
%                   air-vehicle initial conditions, and calls "execute_mission"
%   post_pro_file-  this file will be executed after the mission is run,
%                   used for plotting, scoring, etc.
% 
% OUTPUTS:
%   't_hist' and 'x_hist' timehistory variables are created in the
%   base-workspace.
% 
%   if no post-processing file is provided, "plot_timehist" is assumed
%   and aircraft states are plotted. For no post-processing, use 'none'
%
% EXAMPLE:
%   run_km_sim('p_hyperion.mat','km_checkcase.m');
% 
% see also execute_mission, plot_timehist, plot_trajectory

missionName = mission_file; % Mischa ruggurgi
%crit_info(missionName);

progBar = waitbar(5,'Running','Name','Status'); %ADDED BY MISCHA

% get files if not provided
waitbar(0.15,progBar,'Getting files');
if nargin<3
    post_pro_file = 'plot_timehist'; % basic time history is the default
end
if nargin<2
    vehicle_file = uigetfile('*.mat','Select a vehicle file');
end
if nargin<1
    mission_file = uigetfile('*.m','Select a mission file');
    disp('To run this vehicle/mission combo again without manual selection, use command:');
    disp(['run_km_sim(''', vehicle_file, ''',''', mission_file, ''');']);
end

% load vehicle
waitbar(0.33,progBar,'Loading Vehicle'); %ADDED BY MISCHA
disp('loading vehicle...');
load(vehicle_file,'plane');
current_plane(plane);
disp('vehicle loaded');

% execute mission
waitbar(0.66,progBar,'Executing Mission'); %ADDED BY MISCHA
disp('executing mission...');
run(mission_file);
disp('mission executed');

% post-processing
waitbar(0.9,progBar,'Processing Results'); %ADDED BY MISCHA
run(post_pro_file);

waitbar(1,progBar,'Complete'); %ADDED BY MISCHA
pause(1)
close(progBar)

end


% TODO:   
%       idiot proofing; add truncate/extrapolate option to interpolations
%       instead of just bombing everytime