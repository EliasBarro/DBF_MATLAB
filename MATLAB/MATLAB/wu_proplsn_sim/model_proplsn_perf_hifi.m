function [prop_performance] = model_proplsn_perf_hifi(max_u, density_fs, opt_flg, savepropdata, vehicle_file)

%% model_dyn_prop_perf
% max_u -			max velocity
% density_fs -		freestream density
%  opt_flg     -    0 - no plots, output simple table
%                   1 - just airspeed plots, output simple table
%                   2 - all plots, output complete performance table
% savepropdata -    option to assign propeller performance table off to a plane object
%                   0 - off
%                   1 - use provided "vehicle_name" field
%                   2 - UI to select save
% vehicle_file -  name of vehicle to save off to

%% ===================== performance vs. airspeed =========================
u_fs = (0:1:max_u)';          % m/s

% sim settings
try
    simm = evalin('base','simm');
catch
end
simm.dt    = 0.01;
simm.use_batt  = false;
simm.use_drain = false;
simm.start_dec = 0.8;
simm.output_dt = 0.25;
simm.init_omega_radps = 1500;
assignin('base','simm' ,simm);
settle_time = 0.5;
thr_const   = 1;  % throttle frac
modelname   = 'standalone_propulsion_model';

%overwrite these params just to be safe that they are correct...
options.LoadExternalInput = 'on';
options.ExternalInput = 'throttle, airspeed_m_s, density_kg_m3';
options.StopTime = 'end_time';

%we're going to build up some vectors so suppress the associated warnings 
%#ok<*AGROW>

% define input timeseries - (data, time)
airspeed_data = 0;
airspeed_time = 0;
output_times = [];
for vel_indx = 0:1:max_u
    airspeed_data = [airspeed_data, vel_indx, vel_indx];
    last_time = airspeed_time(end);
    airspeed_time = [airspeed_time, (last_time+simm.dt), (last_time+settle_time)];
    output_times  = [output_times, airspeed_time(end)]; %grab the output at the end of the step
end
end_time = airspeed_time(end);
assignin('base','end_time' ,end_time);

% define input timeseries - (data, time)
assignin('base','airspeed_m_s' ,timeseries(airspeed_data, airspeed_time));
assignin('base','throttle'     ,timeseries([thr_const, thr_const],[0,end_time]));
assignin('base','density_kg_m3',timeseries([density_fs, density_fs],[0,end_time]));

% sim
simOut = sim(modelname,options);

% outputs
angular_rate_pre   = simOut.yout{1}.Values;
electric_power_pre = simOut.yout{2}.Values;
thrust_pre         = simOut.yout{3}.Values;
mech_power_pre     = simOut.yout{4}.Values;
current_pre        = simOut.yout{5}.Values;
torque_pre         = simOut.yout{6}.Values;
batt_state_pre     = simOut.yout{7}.Values;
batt_voltage_pre   = simOut.yout{8}.Values;

% pre allocate
out_size = length(output_times);
angular_rate = zeros(out_size,1);
electric_pwr = zeros(out_size,1);
thrust       = zeros(out_size,1);
mech_power   = zeros(out_size,1);
current      = zeros(out_size,1);
torque       = zeros(out_size,1);
batt_state   = zeros(out_size,1);
batt_voltage = zeros(out_size,1);

% get only the outputs in question...
for indx = 1:1:out_size
    temp_ts           = getsampleusingtime(angular_rate_pre,output_times(indx)); 
    angular_rate(indx)= temp_ts.Data;
    temp_ts           = getsampleusingtime(electric_power_pre,output_times(indx)); 
    electric_pwr(indx)= temp_ts.Data;
    temp_ts           = getsampleusingtime(thrust_pre,output_times(indx)); 
    thrust(indx)      = temp_ts.Data;
    temp_ts           = getsampleusingtime(mech_power_pre,output_times(indx)); 
    mech_power(indx)  = temp_ts.Data;
    temp_ts           = getsampleusingtime(current_pre,output_times(indx)); 
    current(indx)     = temp_ts.Data;
    temp_ts           = getsampleusingtime(torque_pre,output_times(indx)); 
    torque(indx)     = temp_ts.Data;
    temp_ts           = getsampleusingtime(batt_state_pre,output_times(indx)); 
    batt_state(indx)  = temp_ts.Data;
    temp_ts           = getsampleusingtime(batt_voltage_pre,output_times(indx)); 
    batt_voltage(indx)  = temp_ts.Data;
end

%% store off propulsion summary table
%
% this is 1D interpolation: thrust = f(u_fs) and current = f(u_fs)
%
%
prop_performance = [u_fs, thrust, current]; %m/s, N, amps
% save off data to plane object
if savepropdata
    % input from database
    if savepropdata==2
        vehicle_file = uigetfile('*.mat','Select a vehicle file');
    end    
    load(vehicle_file); %plane should now be on the workspace
    disp('vehicle loaded');        
        
    %save
    plane.prop_perform = prop_performance;
%     plane.e_cap = evalin('base',proplsn.max_cap_Ah); % Max energy capacity of battery
    plane_str = ['p_',plane.name,'.mat'];
    save(['database/',plane_str], 'plane');       %this will update teh file
    disp('vehicle record updated with latest propulsion table data');
end


%% Further plotting and analysis as necessary
if opt_flg
    % plot thrust
    figure
    
    subplot(3,2,1);
    plot(u_fs, thrust./9.81, 'b');
    grid on;
    ylabel('thrust (kgf)')
    xlabel('forward velocity (m/s)')
    if max(thrust)>1e-6
        ylim([0, max(thrust./9.81)*1.1]);
    end

    % plot current
    subplot(3,2,2);
    plot(u_fs, current, 'm');
    grid on;
    ylabel('current (A)');
    xlabel('forward velocity (m/s)');
    if max(current)>1e-6
        ylim([0, max(current)*1.1]);
    end

    % plot mechanical power
    subplot(3,2,3);
    plot(u_fs, mech_power, 'cyan');
    grid on;
    ylabel('mech power (W)');
    xlabel('forward velocity (m/s)');
    if max(electric_pwr)>1e-6
        ylim([0, max(electric_pwr)*1.1]);
    end

    % plot power
    subplot(3,2,4);
    plot(u_fs, electric_pwr, 'red');
    grid on;
    ylabel('electric power (W)');
    xlabel('forward velocity (m/s)');
    ylim([0, max(electric_pwr)*1.1]);

    % plot torque
    subplot(3,2,5);
    plot(u_fs, torque, 'black');
    grid on;
    ylabel('torque (Nm)');
    xlabel('forward velocity (m/s)');
    if max(torque)>1e-6
        ylim([0, max(torque)*1.1]);
    end

    % plot angular_rate
    subplot(3,2,6);
    plot(u_fs, angular_rate, 'green');
    grid on;
    ylabel('angular rate (rpm)');
    xlabel('forward velocity (m/s)');
    if max(angular_rate)>1e-6
        ylim([0, max(angular_rate)*1.1]);
    end
    
%     % plot battery voltage
%     figure;
%     plot(u_fs,batt_voltage);
%     grid on;
%     ylabel('voltage');
%     xlabel('forward velocity (m/s)');
%     if max(batt_voltage)>1e-6
%         ylim([0, max(batt_voltage)*1.1]);
%     end

    %% ============== complete performance tables ====================
    if opt_flg==2
        u_const = max_u/2;

        % setup again for throttle:
        thr_vec = -500:100:500;
        num_samps = length(thr_vec);
        thrust = zeros(num_samps,1);
        power  = zeros(num_samps,1);
        torque = zeros(num_samps,1);
        angular_rate     = zeros(num_samps,1);
        current          = zeros(num_samps,1);
        efficiency       = zeros(num_samps,1);
        total_efficiency = zeros(num_samps,1);

        % loop through throttles
        for indx = 1:num_samps
            [thrust(indx), power(indx), torque(indx), angular_rate(indx), current(indx), efficiency(indx),~]...
                = model_proplsn_ss_hifi(u_const, thr_vec(indx), density_fs, settle_time);
        end

        % plot thrust
        figure
        title('performance as it varies with throttle (at half the max airspeed)');
        subplot(3,2,1);
        plot(thr_vec, thrust./9.81, 'b');
        grid on;
        ylabel('thrust (kgf)');
        xlabel('throttle (pwm)');
        ylim([0, max(max(thrust./9.81)*1.1,1e-9)]);

        % plot current
        subplot(3,2,2);
        plot(thr_vec, current, 'm');
        grid on;
        ylabel('current (A)');
        xlabel('throttle (pwm)');
        ylim([0, max(max(current)*1.1,1e-9)]);

        % plot efficiency
        subplot(3,2,3);
        plot(thr_vec, total_efficiency, 'cyan');
        grid on;
        ylabel('total effeciency  (thrust*u/elecpwr)');
        xlabel('throttle (pwm)');
        ylim([0, 1]);

        % plot power
        subplot(3,2,4);
        plot(thr_vec, power, 'red');
        grid on;
        ylabel('mech power (W)');
        xlabel('throttle (pwm)');
        ylim([0, max(max(power)*1.1,1e-9)]);

        % plot torque
        subplot(3,2,5);
        plot(thr_vec, torque, 'black');
        grid on;
        ylabel('torque (Nm)');
        xlabel('throttle (pwm)');
        ylim([0, max(max(torque)*1.1,1e-9)]);

        % plot angular_rate
        subplot(3,2,6);
        plot(thr_vec, angular_rate, 'green');
        grid on;
        ylabel('angular rate (rpm)');
        xlabel('throttle (pwm)');
        ylim([0, max(max(angular_rate)*1.1,1e-9)]);

        %==================== heatmap for total efficiency ======================

        % 2D design space, voltage (x) and airspeed (y)
        n_samps_x = length(thr_vec); %borrow the same vector from before
        n_samps_y = length(u_fs);    %borrow the same vector from before before
        thrust = zeros(n_samps_x, n_samps_y);
        total_efficiency = zeros(n_samps_x, n_samps_y);

        % loop through throttle
        for indx = 1:n_samps_x
            %loop through airspeeds
            for indy = 1:n_samps_y
                %get performance
                [thrust(indx,indy), ~, ~, ~, ~, ~,total_efficiency(indx,indy)]...
                    = model_proplsn_ss_hifi(u_fs(indy), thr_vec(indx), density_fs, settle_time);
            end
		end
		
        % 3D power curve
		if exist('thrust_required.mat','dir')
			thrust_required = [];
			load('thrust_required.mat');     % fills it with the latest "required thrust" for some configuration
			u_tc = thrust_required(:,1);
			thrust_simp = thrust_required(:,2)';
			thrust_req = thrust_simp;
			for indx=2:length(thr_vec)
				thrust_req = [thrust_req; thrust_simp;];
			end
		
			% plot surface of thrust over design space:
			figure
			hold on;
			surf(u_fs, thr_vec, thrust./9.81);
			surf(u_tc, thr_vec(length(thr_vec)/2:length(thr_vec)), thrust_req(length(thr_vec)/2:length(thr_vec),:)./9.81);
			hold off;
			xlabel('velocity of freestream (m/s)');
			ylabel('throttle (pwm)');
			zlabel('thrust (kgf)');
			view(20,20);
			axis vis3d
		end

        % plot surface of efficiency over design space:
        figure
        surf(u_fs, thr_vec, total_efficiency);
        xlabel('velocity of freestream (m/s)');
        ylabel('throttle (pwm)');
        zlabel('total efficiency');
        axis vis3d
        
        %% store off complete propulsion summary table
        %
        % this is 3D interpolation: thrust = f(u_fs, thr, h) and current = f(u_fs, thr, h)
        % prop_performance is a cell array with 5 elements:
        %   - brkpts_u_fs
        %   - brkpts_thr
        %   - brkpts_h
        %   - out_thrust
        %   - out_power
        %
        prop_performance = [u_fs, thrust, current]; %m/s, N, amps
        % save off data to plane object
        if savepropdata
            % input from database
            if savepropdata==2
                vehicle_file = uigetfile('*.mat','Select a vehicle file');
            end    
            load(vehicle_file); %plane should now be on the workspace
            disp('vehicle loaded');
            
            %save
            plane.prop_perform = prop_performance;
        %     plane.e_cap = evalin('base',proplsn.max_cap_Ah); % Max energy capacity of battery
            plane_str = ['p_',plane.name,'.mat'];
            save(['database/',plane_str], 'plane');       %this will update teh file
            disp('vehicle record updated with latest propulsion table data');
        end
    end %end if opt_flg==2
end %end if opt_flg


end % end func
