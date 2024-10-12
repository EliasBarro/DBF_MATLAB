function [density_kg_m3] = model_atmo(h_m)
%MODEL_ATMO determines the air density at this altitude (MSL - meters)
%
% SYNTAX: [density_kg_m3] = model_atmo(h_m)

persistent inda;
if isempty(inda),inda=1;end

air_density_tlu  = [1.225, 1.1673, 1.1117, 1.0581, 1.0065, 0.9569,...
    .9091, .8632, .8191, .7768, .7361, .6971, .6597, .6238, .5895,...
    .5566, .5252, .4951, .4663, .4389, .4127, .3877, .3639, .3108,...
    .2655, .2268, .1937, .1654, .1413, .1207, .1031, .0880];
altitude_bkpts_m = [0:500:11000, 12000:1000:20000];

[density_kg_m3,inda] = lin_interp1(altitude_bkpts_m,air_density_tlu,h_m,inda);

end

