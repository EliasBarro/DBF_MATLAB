function [prop_perf] = get_proplsn_table(all_data, show_data_ext)

%% model_dyn_prop_perf

% *** IT IS EXPECTED THAT AN INIT SCRIPT IS RUN PRIOR TO RUNNING THIS! ****

% we will store the data completely 3 dimensionally
%
% matlab lists 3D data by showing the 2D matrix corresponding to each 3rd
% dimension element. That is, data(:,:,1) is listed first, so they are
% grouped by the last dimension. We will pick this to be density such that
% a table of u & thr can be easily viewed for a given density, this will be
% the slowest varying independent var.
%
% keeping with that theme, flight condition (u_fs) will be varied slower 
% than vehicle setting (thr) so u_fs is the second dimension and thr is the
% third and last.
%
% test points must be conducted sequentially so the 3D table must be
% flattened out to a single vector, ran, and the output data collapsed back
% to a 3D table. 
%
% again, plan is to sweep thru throttle, then move to the next velocity and
% sweep again. After velocity is sweept, density is incremented, until all
% test points are run


% INPUTS:
%   bp_thr      -   breakpoints for throttle, should be between 0 and 1
%   bp_u_fs     -   breakpoints for freestream velocity, in m/s
%   bp_density  -   breakpoints for density, in kg/m^3
%   settle_times-   vector of the settle times between each test point for
%                   a change in each dimension (in order: 1,2,3)
%   all_data    -   if true, will give 4 extra elements in prop_perf
% *** IT IS EXPECTED THAT AN INIT SCRIPT IS RUN PRIOR TO RUNNING THIS! ****

% OUTPUTS:
%   prop_perf   -   cell array containing the breakpoints {1-3} and two 
%                   3D tables of data, {4} containing the thrust at these 
%                   breakpoints, and {5} containing the current in mAh per second

% SYNTAXES:
%   prop_perf = get_proplsn_table(bp_thr, bp_u_fs, bp_density, settle_times, all_data)
%   prop_perf = get_proplsn_table(bp_thr, bp_u_fs, bp_density, settle_times)
%   prop_perf = get_proplsn_table(bp_thr, bp_u_fs, bp_density)
%   prop_perf = get_proplsn_table()


%% SETUP

% start timer
disp('starting propulsion performance analysis...');
tstart = tic;

% defaults
if nargin<2
    show_data_ext = false;
end
if nargin<1
    all_data = false;
end

% get params from proplsn struct
if evalin('base','exist(''proplsn'',''var'')'), proplsn = evalin('base','proplsn');end
if evalin('base','exist(''simm'',''var'')'),    simm    = evalin('base','simm');end
settle_times    = proplsn.settle_times;
bp_dens         = proplsn.bp_dens;
bp_u_fs         = proplsn.bp_u_fs;
bp_thr          = proplsn.bp_thr;

% sim settings
if isfield(proplsn,'dt')
    simm.dt     = proplsn.dt;
else
    simm.dt     = 0.01;
end
simm.use_batt   = true;
simm.use_drain  = false;
simm.start_dec  = 1;
simm.output_dt  = simm.dt;
simm.init_omega_radps = 1500;
assignin('base','simm' ,simm);
modelname       = 'standalone_propulsion_model';


%% INPUT PROCESSING

% DOCUMENTATION - knowing where you are in the "test_points" table.
% Example with nonsensical data, pts_thr = 5, pts_u_fs = 4, pts_dens = >2:
%     (NOTE - this table is transposed to fit the comments easier!)
%  [indx]: 1 2 3 4 5 6...   ...11...  ...16...  ...21...
%  [time]:                  time not important
%  [thr ]: 0 1 2 3 4 0 1 2 3 4 0 1 2 3 4 0 1 2 3 4 0 1 2 3 4 0 1 2 3 4 ...
%  [u_fs]: 0 0 0 0 0 1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 0 0 0 0 0 1 1 1 1 1 ...
%  [dens]: 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 ...
%
%  "minor block" changes occur when mod(indx,minr_blk_sz)==1 [1,6,11,16,21...]
%  "major block" changes occur when mod(indx,majr_blk_sz)==1 [1,21,41,...]
%
%  Note the major block transitions always occur on minor block transitions
%  so the settling time for major blocks (dimension 3) supercedes the
%  settling time for minor blocks (dimension 2), which also supercedes the
%  standard settling time that always happens(dimension 1).
%
%  again, the real table is [time, thr, u_fs, dens] with columns going down

% ----- create the test point table -----
pts_thr     = length(bp_thr);
pts_u_fs    = length(bp_u_fs);
pts_dens    = length(bp_dens);
num_tests   = pts_dens * pts_u_fs * pts_thr;
minr_blk_sz = pts_thr;
majr_blk_sz = pts_thr * pts_u_fs;
test_points = zeros(2*num_tests,4); % time, thr, u, density

% initial test point
test_points(1,1) = 0;
test_points(1,2) = bp_thr(1);
test_points(1,3) = bp_u_fs(1);
test_points(1,4) = bp_dens(1);

% populate rest of test point table
for indx=2:num_tests
    % find dwell time, start with default (which is dim 1)
    dwell_time = settle_times(1);
    
    % check for minor block transition into the last test point
    if mod(indx-1,minr_blk_sz) == 1
        % overwrite dwell time if needed
        dwell_time = settle_times(2);
    end
    
    % check for major block transition into the last test point
    if mod(indx-1,majr_blk_sz) == 1
        % overwrite dwell time if needed
        dwell_time = settle_times(3);
    end
    
    % get location within major and minor blocks
    minor_loc = mod(indx-1,minr_blk_sz)+1;
    major_loc = mod(indx-1,majr_blk_sz)+1;
    
    % define this test point
    test_points(indx,1) = test_points(indx-1,1) + dwell_time;   % time
    test_points(indx,2) = bp_thr(minor_loc);                    % thr
    test_points(indx,3) = bp_u_fs(floor((major_loc-1)/minr_blk_sz)+1);% u_fs
    test_points(indx,4) = bp_dens(floor((indx-1)/majr_blk_sz)+1);     % density
end

% ----- stepify the test_points table -----
% basically, if we make timeseries inputs based on this, we're going to get
% smooth interpolated inputs instead of dwelling at a constant test point.
% We need to doubl the number of points so that things step instead of
% ramp, so where do we step?
%
% well, we want to extract the performance data right before the next test
% point, so we should duplicate the test point data at time (next test
% point - dt)
%
% start at the end so as not to overwrite good data; work backwards. we do
% have to cover the first couple points at the start tho to make sure they
% don't start writing where they need to read

% initialize and cache off starting data
output_times = zeros(num_tests,1);
cache = test_points(1:4,:);

% stepify the final test point first to handle array out of bounds probs
test_points(2*num_tests,:)  = test_points(num_tests,:);
test_points(2*num_tests-1,:)= test_points(num_tests,:);
test_points(2*num_tests,1)  = test_points(num_tests,1)+settle_times(1);
output_times(num_tests)     = test_points(num_tests*2,1)-simm.dt;

% perform stepify operation until you get near indx==1
indx = num_tests-1;
while indx >= 4 
    % copy test data at n out to 2n and 2n-1
    test_points(indx*2,:)   = test_points(indx,:);
    test_points(indx*2-1,:) = test_points(indx,:);
    
    % modify the time at the end of the step to be just prior to the time of the next step
    test_points(indx*2,1)   = test_points(indx+1,1)-simm.dt;
    output_times(indx)      = test_points(indx*2,1);
    
    % decrement indx
    indx = indx - 1;
end

% handle the first 3 data points according to the unaffected cache
for indx=1:3 
    % copy test data at n out to 2n and 2n-1
    test_points(indx*2,:)   = cache(indx,:);
    test_points(indx*2-1,:) = cache(indx,:);
    
    % modify the time at the end of the step to be just prior to the time of the next step
    test_points(indx*2,1)   = cache(indx+1,1)-simm.dt;
    output_times(indx)      = test_points(indx*2,1);
end

% ----- create the timeseries inputs -----
assignin('base','airspeed_m_s' ,timeseries(test_points(:,3),test_points(:,1)));
assignin('base','throttle'     ,timeseries(test_points(:,2),test_points(:,1)));
assignin('base','density_kg_m3',timeseries(test_points(:,4),test_points(:,1)));
end_time = test_points(end,1);
assignin('base','end_time' ,end_time);


%% RUN SIM

%overwrite these params just to be safe that they are correct...
options.LoadExternalInput = 'on';
options.ExternalInput = 'throttle, airspeed_m_s, density_kg_m3';
options.StopTime = 'end_time';

% sim
simOut = sim(modelname,options);

%% OUTPUT PROCESSING

% ----- extract outputs in huge table form -----
angular_rate_pre   = simOut.yout{1}.Values;
electric_pwr_pre   = simOut.yout{2}.Values;
thrust_pre         = simOut.yout{3}.Values;
mech_power_pre     = simOut.yout{4}.Values;
current_pre        = simOut.yout{5}.Values;
torque_pre         = simOut.yout{6}.Values;

% extract only the desired outputs
angular_rate    = interp1(angular_rate_pre.Time , angular_rate_pre.Data , output_times);
electric_pwr    = interp1(electric_pwr_pre.Time , electric_pwr_pre.Data , output_times);
thrust          = interp1(thrust_pre.Time       , thrust_pre.Data       , output_times);
mech_power      = interp1(mech_power_pre.Time   , mech_power_pre.Data   , output_times);
current         = interp1(current_pre.Time      , current_pre.Data      , output_times);
torque          = interp1(torque_pre.Time       , torque_pre.Data       , output_times);

if show_data_ext
    figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(3,1,1);
    plot(thrust_pre.Time,thrust_pre.Data);
    hold on;
    scatter(output_times,thrust);
    legend('thrust','extracted thrust points');
    
    subplot(3,1,2);
    plot(current_pre.Time,current_pre.Data);
    hold on;
    scatter(output_times,current);
    legend('current','extracted current points');
    
    subplot(3,1,3);
    plot(angular_rate_pre.Time,angular_rate_pre.Data);
    hold on;
    scatter(output_times,angular_rate);
    legend('rpm','extracted rpm points');
end

% clear now unneeded copies of data before creating another data copy
clear simOut;
clear angular_rate_pre;
clear electric_power_pre;
clear thrust_pre;
clear mech_power_pre;
clear current_pre;
clear torque_pre;
clear batt_state_pre;

% init 
big0 = zeros(pts_thr,pts_u_fs,pts_dens);
angular_rate_3D = big0;
electric_pwr_3D = big0;
thrust_3D       = big0;
mech_power_3D   = big0;
current_3D      = big0;
torque_3D       = big0;

% ----- construct 3D tables -----
for indx=1:num_tests
    % get dimension location 
    loc_dim1    = mod(indx-1,minr_blk_sz)+1;
    loc_dim2    = floor(((mod(indx-1,majr_blk_sz)+1)-1)/minr_blk_sz)+1; %   o_0   what have I done?! I have created a monster...
    loc_dim3    = floor((indx-1)/majr_blk_sz)+1;
    
    % define this test point
	angular_rate_3D(loc_dim1, loc_dim2, loc_dim3) = angular_rate(indx);
	electric_pwr_3D(loc_dim1, loc_dim2, loc_dim3) = electric_pwr(indx);
	thrust_3D      (loc_dim1, loc_dim2, loc_dim3) = thrust      (indx);
	mech_power_3D  (loc_dim1, loc_dim2, loc_dim3) = mech_power  (indx);
	current_3D     (loc_dim1, loc_dim2, loc_dim3) = current     (indx);
	torque_3D      (loc_dim1, loc_dim2, loc_dim3) = torque      (indx);
end

% if all_data is not set to true, then just output thrust and current
if all_data
    prop_perf = {   bp_thr, bp_u_fs, bp_dens,...
                    angular_rate_3D, electric_pwr_3D, thrust_3D, ...
                    mech_power_3D, current_3D, torque_3D};
else
    eDot_3D   = -current_3D.*1000/3600; % need to convert from A to mAh per sec
    prop_perf = {   bp_thr, bp_u_fs, bp_dens,...
                    thrust_3D, eDot_3D};
end

telapsed = toc(tstart);
disp(['propulsion performance analysis complete in ', num2str(telapsed), 'sec']);

end % end func
