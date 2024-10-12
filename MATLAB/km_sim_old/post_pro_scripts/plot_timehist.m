function plot_timehist()
t       = evalin('base','t_hist');
x_hist  = evalin('base','x_hist');

% x_hist data are [x, v, h, aoa, gamma, psi, e, thr, delta1, delta2, triggr1, triggr2, north, east, mission_phase, mission_phase_step]

%% distances
north   = x_hist(:,13); % meters
east    = x_hist(:,14); % meters
v       = x_hist(:,2);  % meters/sec
h       = x_hist(:,3);  % meters of altitude 
h       = h-h(1);       % change absolute altitude to a clearance altitude

% plot north and east components of groundtrack distance
figure('units','normalized','outerposition',[0 0 1 1]);
ax(1) = subplot(4,2,1);
hold on;
plot(t,north,'blue');grid;
plot(t,east,'green');
legend('North','East');
ylabel('groundtrack distance (m)');

% plot v
ax(2) = subplot(4,2,2);
plot(t,v,'blue');grid;
ylabel('velocity (m/s)');

% plot h
ax(3) = subplot(4,2,3);
plot(t,h,'blue');grid;
ylabel('clearance altitude(m)');


%% angles
aoa     = x_hist(:,4)*180/pi;
gamma   = x_hist(:,5)*180/pi;
psi     = mod(x_hist(:,6),2*pi)*180/pi;

% plot aoa
ax(4) = subplot(4,2,4); hold on;
plot(t,aoa,'blue');grid;
plot(t,gamma,'green');
ylabel('Longitudinal Angles (deg)');
legend('Angle of Attack','Flight Path Angle','Location','Best');

% plot heading
ax(5) = subplot(4,2,5);
plot(t,psi,'blue');grid;
ylabel('Heading (deg)');


%% power system
e       = x_hist(:,7);  % charge state
thr     = x_hist(:,8);  % 0 - 1 fractional thr state

% charge wrt time
ax(6) = subplot(4,2,6);
plot(t,e,'blue');grid;
ylabel('Internal Energy [e]');
ylim([0,max(e)*1.1]);

% throttle wrt time
ax(8) = subplot(4,2,8);
plot(t,thr,'blue');grid;
ylim([0,1]);
ylabel('throttle setting (decimal fraction)');
xlabel('time (sec)');

%% mission phase
ax(7) = subplot(4,2,7);
m_phz = x_hist(:,15);
plot(t,m_phz);grid;
ylabel('mission phase');
xlabel('time (sec)');
yticks([1 2 3 4 5 6 7 8 9 10 11]);
yticklabels({'ground','rotate','transient','cruise-thr','hold-speed','climb-max','glide','opt-climb','opt-glide','turn-gload','turn-radi'});


%% link x axes such that zooming in is coordinated
linkaxes(ax,'x');


end