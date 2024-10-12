function rho = dens(ALT)

%Calulate air density in kg/m^3; see: https://aviation.stackexchange.com/questions/87978/equations-for-computing-temperature-and-air-density-ratios-at-different-altitude

rho = 1.225*(1-(0.0065/288.15)*ALT)^4.25587971;
