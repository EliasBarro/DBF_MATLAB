% Daniel Ruskin - 2023
%% Assumptions/Estimates/Inputs

close all hidden

variables = readmatrix('excel.csv');
% Masses
Struc_mass_min = variables(1); % kg; Minimum structural mass
motor_mass = variables(2); % kg; Estimated mass of the propulsion package

% Course Details
L_ft = variables(3); % ft; Lap length in feet; approximated from the diagram in the rules
L = L_ft/3.281; % m; Lap length in meteres
ALT = variables(4); % m; Course altitude (669 m above SL for Tucson International Modelplex Park Association, TI MPA)
M2_flying_time = variables(5); % seconds; Amount of time for flying M2 (used to calculate # of laps)

% Aircraft Parameters
Cd_min = variables(6); % Estimated minimum drag (AKA parasitic drag); see Snorri Chapter 15 for details
P = variables(7); % ???? Estimated electrical power available at cruise velocity
Cd_i = variables(8); %Induced drage coefficient (taken from GTech 2016 plane)
Cl_to = variables(9); %lift at 0 AoA (taken from GTech 2016 plane)
Cl_max = variables(10); %Max lift at optimal AoA
p = 1.225; %Air density
mu = variables(11); %Ground roll friction coefficient
e_eff = variables(12); %Electrical Efficiency (P_mech/P_elec)
p_eff = variables(13); %Propulsive Efficiency (P_pro/P_mech)
p_eff_to = variables(14); %Propulsive Efficiency at takeoff velocities (P_pro/P_mech) -> propulsive efficeincy decreases with foward velocity (See Snorri pg. 619)
corr = 1; %Correction factor to make equation consistent with 2021
tail_volume = variables(15);

% Constants
g = 9.81; % m/s^2; Gravitational acceleration constant

% Model Varaibles
min_m2_mass = 0.7;%(Struc_mass_min + Prop_mass)*0.3; % kg; Mininum m2 payload
max_m2_mass = 2; % kg; Maximum m2 payload (25 lbs)
n1 = 50; % Number of points to evaluate

min_S = 0.2; % m^2; Minimum wing area
max_S = 0.4; % m^2; Maximum wing area
n2 = 50;

min_AR = 4; % Minimum aspect ratio
max_AR = 7; % Maximum aspect ratio
n3 = 50;

min_ant = 0; % m; Minimum antenna length in M3
max_ant = 1; % m; Maximum antenna length in M3 (1.5 m or ~ 59"), accounting for 62" box
n4 = 50;

% min_power = 500;
% max_power = 1500;
% n5 = 20;
%% Sweep

Sweep_m2_mass = linspace(min_m2_mass, max_m2_mass, n1); % Array of m2 mass values to sweep through
Sweep_S = linspace(min_S, max_S, n2); % Array of wing area values to sweep through
Sweep_AR = linspace(min_AR, max_AR, n3); % Array of AR to sweep through
Sweep_ant = linspace(min_ant, max_ant, n4); % Array of antenna length to sweep through

%% Calculate All-up Weights
%EM = Struc_mass_min + Prop_mass; % kg; Empty mass (i.e. mass without payload)
% prop_mass = 0.1785*Sweep_pow/1000+0.0262;
% EM = Struc_mass_min + repmat(Sweep_S,n2,1) + repmat(prop_mass',1,n5);
EM = Struc_mass_min + motor_mass + Sweep_S;
EW = EM*g;
m2_W = Sweep_m2_mass*g;

%% Calculate M2 Cruise Velocities
t_st_loop=tic;
disp('starting');
% b = ProgressBar([],'IsParallel',true,'Title','running');
% b.setup([],[],[]);

V_max_m2 = zeros(n1,n2,n3); % m/s; Max M2 cruise velocities
counter = zeros(n1,n2,n3); % number of steps needed to converge (for debugging)

parfor i = 1:n1 % Sweep through mass
    
    for j = 1:n2 % Sweep through S
        
        for k = 1:n3 % Sweep through AR

                [V_max_m2(i,j,k), counter(i,j,k)] = perf_V_max(ALT,Sweep_S(j),Cd_min,Sweep_AR(k),EM(j) + Sweep_m2_mass(i),P,e_eff*p_eff); % Calculate max M2 cruise velocities using iterative approach
        end
    end
end

%% Calculate mass of antenna in M3

m_antenna = Sweep_ant*0.2381; % kg; .2381 is linear density of PVC pipe
m3_mass = m_antenna*2; % kg; Added M3 mass is based on 2*(m_antenna) because of counterweight

%% Calculate M3 Cruise Velocities

V_max_m3 = zeros(n4,n1,n2,n3); % m/s; Max M2 cruise velocities
% Switch up matrix dimensions for easier plotting
counter = zeros(n4,n1,n2,n3); % number of steps needed to converge (for debugging)

parfor i = 1:n4 % Sweep through mass
    
    for j = 1:n1 % Sweep through S
        
        for k = 1:n2 % Sweep through AR
            
            for l = 1:n3
                
                    [V_max_m3(i,j,k,l), counter(i,j,k,l)] = perf_V_max(ALT,Sweep_S(k),Cd_min + 1.17*0.021*Sweep_ant(i)/Sweep_S(k),Sweep_AR(l),EM(j) + m3_mass(l),P,e_eff*p_eff); % Calculate max M2 cruise velocities using iterative approach
%                 if V_max_m3(i,j,k,l) > 30
%                     V_max_m3(i,j,k,l) = 30;
%                 end
            end
        end
    end
%     updateParallel();
end
% b.release();

%% Calculate M2 Wing Loading
WL = zeros(n1,n2); % kg/m^2; Wing loading (rows based on m2 payload mass, columns based on S)

parfor i = 1:n1
    for j = 1:n2
        WL(i,j) = (EM(j)+Sweep_m2_mass(i))/Sweep_S(j);
    end
end

%% Calculate Plane Geometry

b = zeros(n3,n2); % m; Wing-span (rows based on AR, columns based on S)
c_ref = zeros(n3,n2);
l_HT = zeros(n3,n2);
S_HT = zeros(n3,n2);
b_HT = zeros(n3,n2);
S_VT = zeros(n3,n2);
b_VT = zeros(n3,n2);


parfor i = 1:n3
    for j = 1:n2
        b(i,j) = sqrt(Sweep_AR(i)*Sweep_S(j));
        c_ref(i,j) = b(i,j)/Sweep_AR(i);
        l_HT(i,j) = sqrt((2*tail_volume*Sweep_S(j)*c_ref(i,j))/(pi*(0.01+0.075)));
        S_HT(i,j) = 0.5*Sweep_S(j)*c_ref(i,j)/l_HT(i,j);
        b_HT(i,j) = sqrt(Sweep_AR(i)/2*S_HT(i,j));
        S_VT(i,j) = 0.04*Sweep_S(j)*b(i,j)/l_HT(i,j);
        b_VT(i,j) = sqrt(2*S_VT(i,j));
    end
end

%% Calculate M1 Score

M1 = zeros(n4,n1,n2,n3);
M1(:) = 1;

%% Calculate M2 Score

M2_raw = zeros(n4,n1,n2,n3);
M2_laps = zeros(n4,n1,n2,n3);

for i = 1:n1
    
    for j = 1:n2
        
        for k = 1:n3

            M2_laps(:,i,j,k) = M2_flying_time/(L/V_max_m2(i,j,k));
            M2_raw(:,i,j,k) = m2_W(i)*M2_laps(1,i,j,k); % Calculate m2 payload weight * # of laps
        
        end
    end
end

%% Calculate M3 Score

M3_raw = zeros(n4,n1,n2,n3);

parfor i = 1:n4
    
    for j = 1:n1
        
        for k = 1:n2
            
            for l = 1:n3
                
            M3_raw(i,j,k,l) = Sweep_ant(i)/(3*L/V_max_m3(i,j,k,l));
            % Calculate antenna length / 3 laps time
            
            end
        end
    end
end


%% Calculate GM Score

GM_raw = zeros(n4,n1,n2,n3);

r = 9.101*10^8; % pa;
r1 = 0.00716; % m;
r2 = 0.00635; % m;

F = zeros(n3,n2); % N; Total test weight
parfor i = 1:n3
    for j = 1:n2
        testweight = abs((2*pi*r*((r2)^4 - (r1)^4))/(b(i,j)*r2));
        if testweight > 600
            testweight = 600;
        end
        F(i,j) = testweight;
    end
end

parfor i = 1:n1
    
    for j = 1:n2
        
        for k = 1:n3
            
            GM_raw(:,i,j,k) = (F(k,j) + EW(j) + Sweep_m2_mass(i))/(EW(j) + Sweep_m2_mass(i)); % Calculate total test weight / max aircraft weight (based on M2 weight)
        
        end
    end
end

allM2 = M2_raw;
allM3 = M3_raw;
allGM = GM_raw;

%% Apply Score Constraint
%Apply box constraint

parfor i = 1:n3
    for j = 1:n2
        if l_HT(i,j) + b_HT(i,j) + b_VT(i,j) > 1.1938       
            M2_raw(:,:,j,i) = NaN;
            M3_raw(:,:,j,i) = NaN;
            GM_raw(:,:,j,i) = NaN;
        end
    end
end

% Apply constraint based on takeoff run for m2 and m3
to_run2 = zeros(n4,n1,n2);
to_run3 = zeros(n4,n1,n2);
parfor i = 1:n4
    for j = 1:n1
        for k = 1:n2
            v_lof2 = 1.556*sqrt((EM(k)+Sweep_m2_mass(j))*9.81/(p*Sweep_S(k)*Cl_max));
            v_lof3 = 1.556*sqrt((EM(k) + m3_mass(i))*9.81/(p*Sweep_S(k)*Cl_max));
            L2 = 0.5*p*Sweep_S(k)*Cl_to*(v_lof2^2)/2;
            L3 = 0.5*p*Sweep_S(k)*Cl_to*(v_lof3^2)/2;
            D2 = 0.5*p*Sweep_S(k)*(Cd_i + Cd_min)*(v_lof2^2)/2;
            D3 = 0.5*p*Sweep_S(k)*(Cd_i + Cd_min + 0.007255060729*Sweep_ant(i))*(v_lof3^2)/2;
            T_eff2 = P*e_eff*p_eff_to*sqrt(2)/v_lof2;
            T_eff3 = P*e_eff*p_eff_to*sqrt(2)/v_lof3;
            a2 = (1/(EM(k)+Sweep_m2_mass(j)))*(T_eff2 - D2 - mu*((EM(k)+Sweep_m2_mass(j))*9.81 - L2));
            a3 = (1/(EM(k) + m3_mass(i)))*(T_eff3 - D3 - mu*((EM(k) + m3_mass(i))*9.81 - L3));

            to_run2(i,j,k) = v_lof2^2/(2*a2)*corr;
            to_run3(i,j,k) = v_lof3^2/(2*a3)*corr;
            if to_run2(i,j,k) > 13.716 || to_run3(i,j,k) > 13.716
                M2_raw(i,j,k,:) = NaN;
                M3_raw(i,j,k,:) = NaN;
                GM_raw(i,j,k,:) = NaN;
            end
        end
    end
end
% M2 weight must be greater than 30% of AUW
parfor i = 1:n1
    for j = 1:n2
        if (Sweep_m2_mass(i)/(EM(j)+Sweep_m2_mass(i))) < 0.3
            M2_raw(:,i,j,:) = NaN;
            M3_raw(:,i,j,:) = NaN;
            GM_raw(:,i,j,:) = NaN;
        end
    end
end

%WL constraint
parfor i = 1:n1
    for j = 1:n2
        if WL(i,j) > 11
            M2_raw(:,i,j,:) = NaN;
            M3_raw(:,i,j,:) = NaN;
            GM_raw(:,i,j,:) = NaN;
        end
    end
end

%% Calculate Total Score

maxM2 = max(M2_raw,[],'All');
maxM3 = max(M3_raw,[],'All');
maxGM = max(GM_raw,[],'All');

M2 = M2_raw/maxM2;
M3 = M3_raw/maxM3;
GM = GM_raw/maxGM;

allM2 = allM2/maxM2;
allM3 = allM3/maxM3;
allGM = allGM/maxGM;

Score = M2 + M3 + GM;
allScore = allM2 + allM3 + allGM;

%% Find indices that produce max score
max_Score = max(Score(:))
Idx = find(Score(:) == max_Score);
[best_ant,best_m2,best_S,best_AR] = ind2sub(size(Score), Idx);
req_cruise = V_max_m3(ind2sub(size(Score), Idx))
req_cruise2 = V_max_m2(best_m2, best_S, best_AR)
best_S_score = Sweep_S(best_S)
best_AR_score = Sweep_AR(best_AR)
best_m2_score = Sweep_m2_mass(best_m2)
best_ant_score = Sweep_ant(best_ant)
tail_arm = l_HT(best_AR, best_S)
tail_span = b_HT(best_AR, best_S)
h_stab_height = b_VT(best_AR, best_S)
best_EW = EM(best_S)
best_WL = WL(best_m2, best_S)
stall = 1.556*sqrt((EM(best_S)+Sweep_m2_mass(best_m2))*9.81/(p*Sweep_S(best_S)*Cl_max))
Thrust = P*e_eff*p_eff_to*sqrt(2)/stall/9.81


%% Sensititivy Plot
fig = figure;
slice_AR = 19;
slice_S = 20;
valid_slice = find(~isnan(Score(:,:,slice_S,slice_AR)));
[i,j] = ind2sub([n4 n1],valid_slice);
valid_ant = Sweep_ant(i)';
valid_m2 = Sweep_m2_mass(j)';
valid_cruise = zeros(length(valid_m2),1);
valid_score = zeros(length(valid_m2),1);
parfor a = 1:length(i)
    valid_cruise(a) = V_max_m3(i(a),j(a),slice_S,slice_AR);
    valid_score(a) = Score(i(a),j(a),slice_S,slice_AR);
end

t_elapsed_loop = toc(t_st_loop);
disp(['complete:, elapsed time: ', num2str(t_elapsed_loop), 'sec']);


xq = unique(Sweep_m2_mass(j));
yq = unique(Sweep_ant(i))';
vq = griddata(valid_m2,valid_ant,valid_cruise,xq,yq);

% chart3d = scatter3(valid_m2, valid_ant, valid_cruise, 20, valid_score, 'filled');
% hold on
chart3d = surf(xq,yq,vq,reshape(valid_score,length(yq),length(xq)));

chart3d.DataTipTemplate.DataTipRows(1).Format = '%.2f'; % x
chart3d.DataTipTemplate.DataTipRows(2).Format = '%.2f'; % y
chart3d.DataTipTemplate.DataTipRows(3).Format = '%.2f'; % z
% chart3d.DataTipTemplate.DataTipRows(4).Format = '%.2f'; % c
chart3d.DataTipTemplate.DataTipRows(1).Label = 'Package Weight'; % x
chart3d.DataTipTemplate.DataTipRows(2).Label = 'Antenna Length'; % y
chart3d.DataTipTemplate.DataTipRows(3).Label = 'M3 Cruise'; % z
% chart3d.DataTipTemplate.DataTipRows(4).Label = 'Score'; % c
colormap parula             
colorbar
ylim([0,1.5])
caxis([2 3])
%set(gca,'XLim',[0.5,2],'YLim',[min_ant, max_ant],'ZLim',[34,40])
ylabel('Antenna Length (m)')
xlabel('M2 Payload (kg)')
zlabel('M3 Cruise (m/s)')
title('Mission Score as a Function of M2 Payload, M3 Antenna Length, and M3 Cruise Velocity')
DragDataTip(fig);

%% Isobar
figure
colormap parula

hold on

c1 = contourf(yq,xq,reshape(valid_score,length(yq),length(xq))',50,'LineColor','none');
c2 = contour(yq,xq,vq','-k',"ShowText",true,"LabelFormat","%0.1f m/s",'LabelSpacing',600);
shading interp
colorbar
clim([1 3])

%% Old Plots
% valid = find(~isnan(Score));
% [i,j,k,l] = ind2sub([n4 n1 n2 n3],valid);
% valid_m2 = Sweep_m2_mass(j)';
% valid_ant = Sweep_ant(i)';
% valid_s = Sweep_S(k)';

% for a = 1:length(i)
%     valid_cruise(a) = V_max_m3(i(a),j(a),k(a),l(a));
%     valid_score(a) = Score(i(a),j(a),k(a),l(a));
% end

% figure
% [X, Y, Z] = meshgrid(Sweep_ant, Sweep_m2_mass, Sweep_S);
% hScatter = scatter3(X(:), Y(:), Z(:), 2, 'filled');
% colormap jet
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
% set(gca,'XLim',[min_ant, max_ant],'YLim',[min_m2_mass, max_m2_mass],'ZLim',[min_S max_S])
% xlabel('Antenna Length (m)')
% ylabel('M2 Payload Mass (kg)')
% zlabel('S (m^2)')
% 
% % set the colorbar
% colorbar
% caxis([floor(0) ceil(3)])

% cycle through the fourth dimension independent variable (i.e., time)
% for k = 10	
%     
%     % visualize the fifth dimension via the marker color
%     C = reshape(Score(:, :, :, k), length(Sweep_ant)*length(Sweep_m2_mass)*length(Sweep_S), 1);  
%        
%     % update the plot
%     set(hScatter, 'CData', C, 'SizeData', 25)  
%     title(['\itScore = \it{f}\rm\bf(\itAntenna Length (m), M2 Payload Mass (kg), S (m^2), AR\rm\bf) \it@ AR = \rm\bf' num2str(Sweep_AR(k))])
%     drawnow
%     
%     % pause for a while
%     pause(3)   
%     
% end