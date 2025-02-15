%Calculate T/W Needed for Takeoff

function TW = takeoff_mod(WL, CLmax, CLto, CDto, SG, mu, ALT)

g = 9.81;

corr = 1; % Correction factor; seems to be needed to make these equations work for RC scale

VLOF = 1.1*sqrt((2*WL)/(dens(ALT)*CLmax));

q = dyn_P(ALT, VLOF/sqrt(2));

TW = corr*(VLOF.^2/(2*g*SG) + q*CDto./WL + mu*(1-q*CLto./WL)); %See Snorri 3-2

end