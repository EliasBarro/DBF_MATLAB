function plot_timehist()
t       = evalin('base','t_hist');
x_hist  = evalin('base','x_hist');
missionName = evalin('caller', 'missionName');

% x_hist data are [x, v, h, aoa, gamma, psi, e, thr, delta1, delta2, triggr1, triggr2, north, east, mission_phase, mission_phase_step]

%% distances
north   = x_hist(:,13); % meters
east    = x_hist(:,14); % meters
v       = x_hist(:,2);  % meters/sec
h       = x_hist(:,3);  % meters of altitude 
h       = h-h(1);       % change absolute altitude to a clearance altitude

% Official WashU Color Pallette
washuGreen = [0.19 0.34 0.21];
washuRed = [0.67 0.15 0.20];

% plot north and east components of groundtrack distance
figure('units','normalized','outerposition',[0 0 1 1]);
sgtitle(missionName); % mischa hrvhv
ax(1) = subplot(4,2,1);
hold on;
axis tight
plot(t,north,Color=washuRed,LineWidth=1.5);grid;
plot(t,east,Color=washuGreen,LineWidth=1.5);
legend('North','East');
ylabel('groundtrack distance (m)');

% plot v
ax(2) = subplot(4,2,2);
plot(t,v,Color = washuRed,LineWidth=1.5);grid;
axis tight
ylabel('velocity (m/s)');

% Velocity Orange = [0.85, 0.325, 0.098]

% plot h
ax(3) = subplot(4,2,3);
plot(t,h,Color=washuRed,LineWidth=1.5);grid;
axis tight
ylabel('clearance altitude(m)');

% Altitude Navy = [0 0 0.39]

%% angles
aoa     = x_hist(:,4)*180/pi;
gamma   = x_hist(:,5)*180/pi;
psi     = mod(x_hist(:,6),2*pi)*180/pi;

% plot aoa
ax(4) = subplot(4,2,4); hold on;
plot(t,aoa,Color = washuRed,LineWidth=1.5);grid;
plot(t,gamma,Color=washuGreen,LineWidth=1.5);
axis tight
ylabel('Longitudinal Angles (deg)');
legend('Angle of Attack','Flight Path Angle','Location','Best');

% plot heading
ax(5) = subplot(4,2,5);
plot(t,psi,Color=washuGreen,LineWidth=1.5);grid;
axis tight
ylabel('Heading (deg)');


%% power system
e       = x_hist(:,7);  % charge state
thr     = x_hist(:,8);  % 0 - 1 fractional thr state

% charge wrt time
ax(6) = subplot(4,2,6);
plot(t,e,Color=washuRed,LineWidth=1.5);grid;
axis tight
ylabel('Internal Energy [e]');
ylim([0,max(e)*1.1]);

% throttle wrt time
ax(8) = subplot(4,2,8);
plot(t,thr,Color=washuGreen,LineWidth=1.5);grid;
axis tight
ylim([0,1]);
ylabel('throttle setting (decimal fraction)');
xlabel('time (sec)');

%% mission phase
ax(7) = subplot(4,2,7);
m_phz = x_hist(:,15);
plot(t,m_phz,Color =washuRed,LineWidth=1.5);grid;
ylabel('mission phase');
axis tight
xlabel('time (sec)');
yticks([1 2 3 4 5 6 7 8 9 10 11]);
yticklabels({'ground','rotate','transient','cruise-thr','hold-speed','climb-max','glide','opt-climb','opt-glide','turn-gload','turn-radi'});


%% link x axes such that zooming in is coordinated
linkaxes(ax,'x');

end