% check for plane obj
if ~exist('plane','var')
    error('no plane object found in current workspace!');
end

% check for common issues
if any(any(isnan(plane.aero)))
    warning('aero table has NaN!!!');
end
if isfield(plane,'alt_aero')
    if any(any(isnan(plane.alt_aero)))
        warning('alt aero table has NaN!!!');
    end
    if any(plane.alt_aero(:,3)<0) || any(plane.aero(:,3)<0)
        warning('negative drag in initial table data!');
    end
else
    if any(plane.aero(:,3)<0)
        warning('negative drag in initial table data!');
    end
end
%% plot CL and CD vs alpha
% CL v alpha
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1);
hold on;
plot(plane.aero(:,1)*180/pi,plane.aero(:,2),'LineWidth',2);
if isfield(plane,'alt_aero')
    plot(plane.alt_aero(:,1)*180/pi,plane.alt_aero(:,2),'LineStyle','--');
    legend('primary','alternate');grid;
else
    plot(plane.aero(:,1)*180/pi,plane.aero(:,2)); grid;
end
xlabel('alpha (deg)');
ylabel('CL');

% CD v alpha
subplot(1,2,2);
hold on;
plot(plane.aero(:,1)*180/pi,plane.aero(:,3),'LineWidth',2);
if isfield(plane,'alt_aero')
    plot(plane.alt_aero(:,1)*180/pi,plane.alt_aero(:,3),'LineStyle','--');
    legend('primary','alternate');grid;
else
    plot(plane.aero(:,1)*180/pi,plane.aero(:,3)); grid;
end
xlabel('alpha (deg)');
ylabel('CD');


%% plot all aero tables CLCD
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,1,1);
hold on;
plot(plane.aero(:,2),plane.aero(:,3),'LineWidth',2);
if isfield(plane,'alt_aero')
    plot(plane.alt_aero(:,2),plane.alt_aero(:,3),'LineStyle','--');
    plot([min(plane.alt_aero(:,2)),max(plane.alt_aero(:,2))],[0,0],'--k');
    legend('primary','alternate');grid;
else
    plot([min(plane.aero(:,2)),max(plane.aero(:,2))],[0,0],'--k'); grid;
end
xlabel('CL');
ylabel('CD');


%% plot thrust table near 2000m
% setup
prop = plane.prop_perf;
indx_dens = length(prop{3})-1;
thrust = prop{4};
thrust = thrust(:,:,indx_dens);
thrust = thrust';

% thrust plot
subplot(2,2,3);
surface(prop{1},prop{2},thrust);
xlabel('throttle decimal');
ylabel('airspeed (m/s)');
zlabel('thrust (N)');grid;
title(['thrust table at density ', num2str(prop{3}(indx_dens)), ' kg/m3']);
view([210 40]);

%% current table
% setup
currnt = prop{5};
currnt = currnt(:,:,indx_dens);
currnt = -currnt';

% current plot
subplot(2,2,4);
surface(prop{1},prop{2},currnt);
xlabel('throttle decimal');
ylabel('airspeed (m/s)');
zlabel('Energy Rate ([e]/s)');grid;
title(['Energy Rate table at density ', num2str(prop{3}(indx_dens)), ' kg/m3']);
view([210 40]);