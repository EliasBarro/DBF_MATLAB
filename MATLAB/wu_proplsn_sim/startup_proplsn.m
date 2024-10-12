% startup script - if you start matlab in this folder, it will find this
% file and run it automatically, setting up this sim to work properly.
%init file

path = fileparts(mfilename('fullpath')); %identifies path to this file

% add paths
%addpath(fullfile(path, '../../../SIM_DEV/proplsn_sim'));
addpath(fullfile(path,'UIUCpropDB'));
addpath(fullfile(path,'proplsn_init_scripts'));
addpath(fullfile(path,'first_order_motor_model'));
addpath(path);
%rmpath(path);


% show init complete
disp('[paths defined] - proplsn sim');
clear path;