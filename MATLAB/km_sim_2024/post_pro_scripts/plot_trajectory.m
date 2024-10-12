function plot_trajectory()
x_hist  = evalin('base','x_hist');

% states are [x, v, h, aoa, gamma, psi, e, thr, delta_tr, trigger_tr, north, east, mission_phase, mission_phase_step];
h       = x_hist(:,3);
north   = x_hist(:,13);
east    = x_hist(:,14);

% convert to clearance altitude
h = h-h(1);

% plot groundtrack
figure('units','normalized','outerposition',[0.25 0 0.5 1]);
plot3(east, north, h,'LineWidth',2,'Color',[10,96,150]/255);
xlabel('East (m)');
ylabel('North(m)');
zlabel('Clearance Altitude (m)');grid;
axis vis3d;
daspect([1 1 1]);
view(0,90);

end