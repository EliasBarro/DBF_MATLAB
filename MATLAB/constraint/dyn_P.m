% Calculate Dynamic Pressure

function q = dyn_P(ALT, V)

% T = (1-0.0000068756*ALT); % Calculate ISA temperature ratio at given altitude; see Snorri Example 16-2

% rho = 1.225*T^4.2561; %Calulate air density in kg/m^3; see: https://aviation.stackexchange.com/questions/87978/equations-for-computing-temperature-and-air-density-ratios-at-different-altitude

rho = 1.225*(1-(0.0065/288.15)*ALT).^4.25587971; %Calulate air density in kg/m^3; see: https://aviation.stackexchange.com/questions/87978/equations-for-computing-temperature-and-air-density-ratios-at-different-altitude

q = 0.5*rho.*V.^2;


