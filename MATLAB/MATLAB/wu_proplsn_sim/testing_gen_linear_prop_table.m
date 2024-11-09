% cleanup
clear; close all; clc;

% call the propulsion init script
sl_static_thrust = 10;
sl_static_edot = 5;
zero_thrust_speed = 30;
max_alt = 500;
prop_perf = gen_linear_prop_table(sl_static_thrust, sl_static_edot, zero_thrust_speed, max_alt);

% extract breakpoints and tables
bp_thr      = prop_perf{1};
bp_u_fs     = prop_perf{2};
bp_dens     = prop_perf{3};
thrust_3D   = prop_perf{4};
current_3D  = prop_perf{5};

% convert thrust from N to kgf
thrust_3D = thrust_3D./9.81;

% define density points
dens_pts = [1,numel(bp_dens)];

% plot a couple thrust tables at specific densities 
for dens_indx = dens_pts
    figure;
    [surf_X, surf_Y] = meshgrid(bp_u_fs, bp_thr);
    surf_Z = thrust_3D(:,:,dens_indx);
    surface(surf_X, surf_Y, surf_Z);
    xlabel('velocity of freestream (m/s)');
    ylabel('throttle (pwm)');
    zlabel('thrust (kgf)');
    title(['Thrust Performance at density = ', num2str(bp_dens(dens_indx)), ' kg/m^3']);
    view(20,20);
    axis vis3d
end

% plot a couple current tables at specific densities 
for dens_indx = dens_pts
    figure;
    [surf_X, surf_Y] = meshgrid(bp_u_fs, bp_thr);
    surf_Z = current_3D(:,:,dens_indx);
    surface(surf_X, surf_Y, surf_Z);
    xlabel('velocity of freestream (m/s)');
    ylabel('throttle (pwm)');
    zlabel('current (A)');
    title(['Current at density = ', num2str(bp_dens(dens_indx)), ' kg/m^3']);
    view(20,20);
    axis vis3d
end