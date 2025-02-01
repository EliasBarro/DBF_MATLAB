function plot3D()
t       = evalin('base','t_hist');
x_hist  = evalin('base','x_hist');

% Official WashU Color Pallette
washuGreen = [0.19 0.34 0.21];
washuRed = [0.67 0.15 0.20];

% x_hist data are [x, v, h, aoa, gamma, psi, e, thr, delta1, delta2, triggr1, triggr2, north, east, mission_phase, mission_phase_step]

%% distances
north   = x_hist(:,13); % meters 
east    = x_hist(:,14); % meters
v       = x_hist(:,2);  % meters/sec
h       = x_hist(:,3);  % meters of altitude
h       = h-h(1);       % change absolute altitude to a clearance altitude

C=[v, v, v];

C(:,1) = v-min(v);
C(:,2) = v-min(v);

C=C/max(C);

C(:,3) = 0;

north = north * 3.281;  % convert from meters to feet
east = east * 3.281;
h = h * 3.281;

for i = 1:28 %iterate through phases
    color = washuGreen;

    %selected pahses
    if i == 24
        color = [1 0 0]; %red
    end
    if i == 25
        color = [0 0 1]; %green
    end
     if i == 26
        color = [0 0 1]; %blue
    end
  
    phase = (x_hist(:,16) == i);
    if ~isempty(phase) 
        scatter3(north(phase), east(phase), h(phase), 5, C(phase,:), MarkerEdgeColor=color);
        axis equal;
        hold on
    end
end


disp("Total flight time: " + max(t));

end