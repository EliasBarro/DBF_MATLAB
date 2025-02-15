close all

%% Plotting Settings
c = [0.647 0.078 0.090;
     0.424 0.451 0.451;
     0.784 0.784 0.784;
     0 0.451 0.376;
     0 0.373 0.522;
     1 0.8 0];
set(0,'DefaultAxesColorOrder',c)
set(0,'defaultLineLineWidth',1.5)

%% Constraints / Estimated Values
%2025 - should I test this at different takeoff values? At the predicted
%takeoff value?
V=31.3; %m/s

AR = 4.1163; %Aspect ratio
AUW = 8.0175; %kg
CDmin = 0.08; %Minimum drag coefficient
CLto = 1; %Coefficient of lift at takeoff
CDto = 0.08; %Coefficient of drag at takeoff
CLmax = 1.5;

mu = 0.00; %Ground roll friction coefficient
e_eff = 0.9; %Electrical Efficiency (P_mech/P_elec)
p_eff = 0.83; %Propulsive Efficiency (P_pro/P_mech)
p_eff_to = 0.35; %Propulsive Efficiency at takeoff velocities (P_pro/P_mech) -> propulsive efficeincy decreases with foward velocity (See Snorri pg. 619)

SG = 7.62; %ft->m
ALT_ground = 668; %m
ALT_cruise = 768; %m

max_TW = 1.75;
max_Power = 1000; %W

min_WL = 4; %kgf/m^2

shade_on = 0; % Turn shading on and off; off is easier for analysis since it allows you to actually click the data points, but shading is pretty for the report

%% Sweep

WL_kg = 2:0.1:15; %kg/m^2
WL = WL_kg*9.81;

%% Calculate T/W

TW_to = takeoff_mod(WL, CLmax, CLto, CDto, SG, mu, ALT_ground);

TW_c = cruise(WL, AR, CDmin, V, ALT_cruise);

TW_tot = TW_to + TW_c;

minTW = find(TW_tot == min(TW_tot));

%% Calculate Power

W_to = reqpower_mod(TW_to, AUW, WL, CLmax, e_eff, p_eff_to, ALT_ground);
W_c = reqpower(TW_c, AUW, V, e_eff, p_eff);

%% Plot

if shade_on ~= 1
    
    figure
    hold on
    plot(WL_kg, TW_to);
    plot(WL_kg, TW_c);

    xlabel('Wing-Loading (kg/m^2)');
    xline(4, "--");
    xline(10, '--');
    yline(2.25, '--');
    ylabel('Required Thrust-to-Weight');
    p = plot(9.9,0.944911,'rp','MarkerSize',15,'MarkerFaceColor','y','MarkerEdgeColor','#EDB120');
    text(10.2,0.9,['(' num2str(9.9) ' kg/m^2, ' num2str(0.944911) ')'])
    hold off
    legend('25 ft Takeoff','31.3 m/s Cruise');
    
    figure
    hold on
    plot(WL_kg, W_to);
    plot(WL_kg, W_c);

        xline(4, "--");
    xline(10, '--');
    yline(2000,'--');
    xlabel('Wing-Loading (kg/m^2)');
    ylabel('Required Power (W)');
    plot(9.9,1948.73,'rp','MarkerSize',15,'MarkerFaceColor','y','MarkerEdgeColor','#EDB120')
    text(10.2,1900,['(' num2str(9.9) ' kg/m^2, ' num2str(1948.73) ' W)'])
    hold off
    ylim([0,2500])
    legend('25 ft Takeoff','31.3 m/s Cruise');
    
else
    
    figure;
    shade(WL_kg,TW_to,WL_kg,TW_c,'LineWidth',1.5,'FillType',[1 0; 2 0],'FillColor','k');
    xlim([min(WL_kg) 15]);
    hold on
    xline(min_WL,'--');
    yline(max_TW,'--');
    xlabel('Wing-Loading [kgf/m^2]');
    ylabel('Required Thrust-to-Weight [-]');
    
    boxtop = ceil(max([max(TW_to,[],'All'),max(TW_c,[],'All')]));
    
    f1 = [1 2 3 4];
    v1 = [0 0; 0 boxtop; min_WL boxtop; min_WL 0];
    
    patch('Faces',f1,'Vertices',v1,'FaceColor','k','FaceAlpha',.2,'EdgeColor','none');
    
    f2 = [1 2 3 4];
    v2 = [0 max_TW; 0 boxtop; max(WL_kg) boxtop; max(WL_kg) max_TW];
    patch('Faces',f2,'Vertices',v2,'FaceColor','k','FaceAlpha',.2,'EdgeColor','none');
    
    legend('25 ft Takeoff','22.5 m/s Cruise');
    
    figure
    shade(WL_kg,W_to,WL_kg,W_c,'LineWidth',2,'FillType',[1 0; 2 0],'FillColor','k');
    xlim([min(WL_kg) 15]);
    hold on
    xline(min_WL,'--');
    yline(max_Power,'--');
    xlabel('Wing-Loading [kgf/m^2]');
    ylabel('Required Power [W]');
    
    boxtop = ceil((max([max(W_to,[],'All'),max(W_c,[],'All')])/1000))*1000;
    
    f1 = [1 2 3 4];
    v1 = [0 0; 0 boxtop; min_WL boxtop; min_WL 0];
    
    patch('Faces',f1,'Vertices',v1,'FaceColor','k','FaceAlpha',.2,'EdgeColor','none');
    
    f2 = [1 2 3 4];
    v2 = [0 max_Power; 0 boxtop; max(WL_kg) boxtop; max(WL_kg) max_Power];
    patch('Faces',f2,'Vertices',v2,'FaceColor','k','FaceAlpha',.2,'EdgeColor','none');
    legend('25 ft Takeoff','22.5 m/s Cruise');
    
end