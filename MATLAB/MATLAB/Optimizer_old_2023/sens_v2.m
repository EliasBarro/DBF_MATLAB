% Michael Qiu - 2022

%% Assumptions/Estimates/Inputs

close all hidden

% Masses
Pack_mass = 0.23; % kg; Mass of a package
Syr_mass = 0.02; % kg; Mass of a syringe
Struc_mass_min = 0.5; % kg; Minimum structural mass
Struc_mass_per_pack = 1; % kg; Estimated additional structural mass required for each additional package (i.e. fuse + wing need to get bigger)
Prop_mass = 1.2; % kg; Estimated mass of the propulsion package;

% Course Details
L_ft = 2500; % ft; Lap length in feet; approximated from the diagram in the rules
L = L_ft/3.281; % m; Lap length in meters
ALT = 400; % m; Course altititude (400 m above SL for Wichita)

% Aircraft Parameters
WL = 6.5; % kg/m^2; Estimated wing loading needed for 25 ft TOFL
Cd_min = 0.15; % Estimated minimum drag (AKA parasitic drag); see Snorri Chapter 15 for details
AR = 5.5; % Estimated Aspect Ratio (Doesn't really change anything velocity wise if kept between 4-7)
P = 1500; % W; Estimated electrical power avaliable at cruise velocity
e_eff = 0.9; % Estimated electrical efficiency (Electrical power -> Mechanical power)
p_eff = 0.78; % Estimated propulsive efficiency (Mechanical power -> Propulsive power)
Max_b = 8/3.281; % m; Max wing-span

% Ground Mission Timing
Time2box = 5; % s; Estimated time to run-to/back-from the mission box during the ground mission
Time2loadsyringe = 10; % s; Estimated time to load 10 syringes
Time2unloadsyringe = 5; % s; Estimated time to unload 10 syringes
Time2loadpackage = 10; % s; Estimated time to load 1 pacakge

% Constants
g = 9.81; % m/s^2; Gravitational acceleration constant

% Model Variables
min_deploy = 1; % Min number of deployments
max_deploy = 7; % Max number of deployments to sweep through
n = max_deploy^2; % Number of points to evaluate

%% Toggles

%% Sweep

N_deploy = linspace(min_deploy,max_deploy,n); % Number of deployments in M3 
N_syringes = 10*N_deploy; % Number of syringes in M2

%% Calculate All-up Weights and Wing Areas

EM = Struc_mass_min + Struc_mass_per_pack*N_deploy + Prop_mass; % kg; Empty mass (i.e. mass without payload);
M = N_syringes*Syr_mass + EM; % kg; All-up mass
W = M*g; % N; All-up weight (AUW)

S = M/WL; % m^2; Wing area for each AUW assuming constant WL (note: since WL is in kg/m^2 we divide all-up mass not AUW)

%% Calculate M2 Cruise Velocities

V_max = zeros(1,n); % m/s; Max M2 cruise velocities
counter = zeros(1,n); % number of steps needed to converge (for debugging)

for i = 1:n

    [V_max(i), counter(i)] = perf_V_max(ALT,S(i),Cd_min,AR,M(i),P,e_eff*p_eff); % Calculate max M2 cruise velocities using iterative approach
    
end

%% Calculate M1 Score

M1 = zeros(1,n);
M1(:) = 1;

%% Calculate M2 Score

M2_raw = zeros(1,n);

for j = 1:n
    
    M2_raw(j) = N_syringes(j)/(3*L/V_max(j));
    
end

M2_norm = M2_raw/max(M2_raw,[],'All');
M2 = 1 + M2_norm;

%% Calculate M3 Score

M3 = 2 + N_deploy/(max(N_deploy,[],'All'));

%% Calculate GM Score

GM_time = (Time2box*4 + (N_syringes(j)/10)*Time2loadsyringe + (N_syringes/10)*Time2unloadsyringe + N_deploy*Time2loadpackage);

GM = min(GM_time,[],'All')./GM_time;

%% Calculate Total Score

Score = M1 + M2 + M3 + GM;

%% Calculate Span

b = sqrt(AR*S); % m; Wing-span based on standard mean chord

%% Plot

nIDs = 4;
alphabet = ('a':'z').';
chars = num2cell(alphabet(1:nIDs));
chars = chars.';
charlbl = strcat('(',chars,')');
posx = 0.025;
posy = 0.9;

figure
subplot(2,3,1)
xlabel('M2 #of Syringes');
xlim([5 75]);
yyaxis right
plot(N_syringes,V_max,'Linewidth', 2);
ylabel('M2 Max Cruise Velocity [m/s]','Color',[17 17 17]/255);
ylim([18 32]);
hold on
yyaxis left
plot(N_syringes,M2,'Linewidth', 2);
ylim([1.2 2.05]);
text(posx,posy,charlbl{1},'Units','normalized','FontSize',16)
hold off
xline(40,'--');
ylabel('M2 Score','Color',[17 17 17]/255);
% ax = gca;
% ax.YAxis(1).Color = '#a51417';
% ax.YAxis(2).Color = '#6c7373';
xticks(10:10:70);

subplot(2,3,2)
plot(N_deploy,M3,'Linewidth', 2);
text(posx,posy,charlbl{2},'Units','normalized','FontSize',16)
xline(4,'--');
ylabel('M3 Score');
xlabel('M3 #of Deployments');
xlim([0.5 7.5]);
ylim([2.1 3.05]);
xticks(1:1:7);

subplot(2,3,3)
plot(N_deploy,GM,'Linewidth', 2);
text(posx,posy,charlbl{3},'Units','normalized','FontSize',16)
xline(4,'--');
ylabel('GM Score');
xlabel('M3 #of Deployments');
xlim([0.5 7.5]);
ylim([0.5 1.04]);
xticks(1:1:7);

subplot(2,3,[4 5 6])
plot(N_deploy,Score,'Linewidth', 2);
text(posx,posy,charlbl{4},'Units','normalized','FontSize',16)
xline(4,'--');
ylabel('Total Mission Score');
xlabel('M3 #of Deployments');
xlim([0.5 7.5]);
ylim([5.2 6.6]);
xticks(1:1:7);

figure
subplot(2,2,1)
plot(N_deploy,M,'Linewidth', 2);
text(posx,posy,charlbl{1},'Units','normalized','FontSize',16)
ylabel('All-up Weight [kgf]');
xlabel('M3 #of Deployments');
xlim([0.5 7.5]);
ylim([2.7 10.3]);
xticks(1:1:7);

subplot(2,2,2)
plot(N_deploy,S,'Linewidth', 2);
text(posx,posy,charlbl{2},'Units','normalized','FontSize',16)
ylabel('Wing Area [m^2]');
xlabel('M3 #of Deployments');
xlim([0.5 7.5]);
ylim([0.4 1.6]);
xticks(1:1:7);

subplot(2,2,3)
plot(N_deploy,GM_time,'Linewidth', 2);
text(posx,posy,charlbl{3},'Units','normalized','FontSize',16)
ylabel('GM Time [s]');
xlabel('M3 #of Deployments');
xlim([0.5 7.5]);
ylim([100 200]);
xticks(1:1:7);

subplot(2,2,4)
plot(N_deploy,b,'Linewidth', 2);
yline(Max_b,'Color','#D95319');
text(posx,posy,charlbl{4},'Units','normalized','FontSize',16)
xline(4,'--');
ylabel('Wing Span [m]');
xlabel('M3 #of Deployments');
xlim([0.5 7.5]);
ylim([1.5 3]);
xticks(1:1:7);

% figure
% plot(N_deploy,Score,N_deploy,M3,N_deploy,M2,N_deploy,M1,N_deploy,GM);