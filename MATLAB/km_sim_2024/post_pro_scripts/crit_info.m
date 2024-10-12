%function crit_info(x_hist)
function crit_info(x_hist,t_hist) %Mischa added t_hist

x=x_hist(1,3);
increment = 1;
laps = 0; % MISCHA
%angle = x_hist(1,6);
n = 1; % MISCHA
maxNorth = max(x_hist(:,13));
minNorth = min(x_hist(:,13));
    
while(x_hist(increment,3)==x)
    increment = increment + 1;
end

% MISCHA:

while(n < length(x_hist(:,13)))
    if (x_hist(n,13) >= maxNorth-0.05)
        laps = laps + 1;
        laps_n = n;
        %finalLapTime = t_hist(n);
        n = n+5;
    else
        n = n+1;
    end
end

n = laps_n;
while(x_hist(n,13) >= minNorth+0.05)
    n = n+1;
end

% MISCHA
finalLapTime = t_hist(n);
totalTime = t_hist(length(t_hist));
time_after_last_lap = totalTime - finalLapTime;

print=["takeoff distance (m): ", x_hist(increment,14)];
disp(print);
print3 = ["takeoff distance (ft)" ,x_hist(increment,14)*3.28084];
disp(print3)

print2=["max cruise speed: ",max(x_hist(:,2))];
disp(print2);

% MISCHA
print4 = ["number of laps: ", laps];
disp(print4);

print5 = ["time at end of last lap (s): ", finalLapTime];
disp(print5);

% print6 = ["time after last lap (s): ", time_after_last_lap];
% disp(print6);