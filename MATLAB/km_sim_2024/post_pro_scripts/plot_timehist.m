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

% plot north and east components of groundtrack distance
figure('units','normalized','outerposition',[0 0 1 1]);
sgtitle(missionName); % mischa hrvhv
ax(1) = subplot(4,2,1);
hold on;
<<<<<<< Updated upstream
axis tight
plot(t,north,Color=[0.67 0 0.14],LineWidth=1.5);grid;
plot(t,east,Color=[0 0.5 0],LineWidth=1.5);
legend('North','East');
ylabel('groundtrack distance (m)');

% WashU Red = [0.67 0 0.14]
% WashU Green = [0 0.5 0]

% plot v
ax(2) = subplot(4,2,2);
plot(t,v,Color = [0.8500, 0.3250, 0.0980],LineWidth=1.5);grid;
axis tight
ylabel('velocity (m/s)');

% Velocity Orange = [0.85, 0.325, 0.098]

% plot h
ax(3) = subplot(4,2,3);
plot(t,h,Color=[0.0 0 0.39]);grid;
axis tight
ylabel('clearance altitude(m)');

% Altitude Navy = [0 0 0.39]
=======
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

>>>>>>> Stashed changes

%% angles
aoa     = x_hist(:,4)*180/pi;
gamma   = x_hist(:,5)*180/pi;
psi     = mod(x_hist(:,6),2*pi)*180/pi;

% plot aoa
ax(4) = subplot(4,2,4); hold on;
<<<<<<< Updated upstream
plot(t,aoa,Color = [0.67 0 0.14]);grid;
plot(t,gamma,Color=[0 0.5 0]);
axis tight
=======
plot(t,aoa,'blue');grid;
plot(t,gamma,'green');
>>>>>>> Stashed changes
ylabel('Longitudinal Angles (deg)');
legend('Angle of Attack','Flight Path Angle','Location','Best');

% plot heading
ax(5) = subplot(4,2,5);
<<<<<<< Updated upstream
plot(t,psi,Color=[0.4940 0.1840 0.5560]);grid;
axis tight
=======
plot(t,psi,'blue');grid;
>>>>>>> Stashed changes
ylabel('Heading (deg)');


%% power system
e       = x_hist(:,7);  % charge state
thr     = x_hist(:,8);  % 0 - 1 fractional thr state

% charge wrt time
ax(6) = subplot(4,2,6);
<<<<<<< Updated upstream
plot(t,e,Color=[0.9290 0.6940 0.1250]);grid;
axis tight
=======
plot(t,e,'blue');grid;
>>>>>>> Stashed changes
ylabel('Internal Energy [e]');
ylim([0,max(e)*1.1]);

% throttle wrt time
ax(8) = subplot(4,2,8);
plot(t,thr,'blue');grid;
<<<<<<< Updated upstream
axis tight
=======
>>>>>>> Stashed changes
ylim([0,1]);
ylabel('throttle setting (decimal fraction)');
xlabel('time (sec)');

%% mission phase
ax(7) = subplot(4,2,7);
m_phz = x_hist(:,15);
<<<<<<< Updated upstream
plot(t,m_phz,Color =[0.74 0.25 0.53]);grid;
ylabel('mission phase');
axis tight
=======
plot(t,m_phz);grid;
ylabel('mission phase');
>>>>>>> Stashed changes
xlabel('time (sec)');
yticks([1 2 3 4 5 6 7 8 9 10 11]);
yticklabels({'ground','rotate','transient','cruise-thr','hold-speed','climb-max','glide','opt-climb','opt-glide','turn-gload','turn-radi'});


%% link x axes such that zooming in is coordinated
linkaxes(ax,'x');

<<<<<<< Updated upstream
=======

>>>>>>> Stashed changes
end