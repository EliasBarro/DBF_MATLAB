loadPlane = 'PhlyzerM1.mat'; 
savePlane = 'PhlyzerM1.mat';
load(loadPlane);

staticThrustNewtons = 106;
static_edot = -27.77;  %Should be negative
pitchSpeed = 26;
maxAlt = 769;

plane.prop_perf = gen_linear_prop_table(staticThrustNewtons, static_edot, pitchSpeed, maxAlt);
save(savePlane, 'staticThrustNewtons', 'static_edot', 'pitchSpeed', 'maxAlt');

%%
loadPlane = 'PhlyzerM2.mat'; 
savePlane = 'PhlyzerM2.mat';
load(loadPlane);

staticThrustNewtons = 106;
static_edot = -27.77;  %Should be negative
pitchSpeed = 26;
maxAlt = 769;

plane.prop_perf = gen_linear_prop_table(staticThrustNewtons, static_edot, pitchSpeed, maxAlt);
save(savePlane, 'staticThrustNewtons', 'static_edot', 'pitchSpeed', 'maxAlt');

%%
loadPlane = 'PhlyzerM3.mat'; 
savePlane = 'PhlyzerM3.mat';
load(loadPlane);

staticThrustNewtons = 106;
static_edot = -27.77;  %Should be negative
pitchSpeed = 26;
maxAlt = 769;

plane.prop_perf = gen_linear_prop_table(staticThrustNewtons, static_edot, pitchSpeed, maxAlt);
save(savePlane, 'staticThrustNewtons', 'static_edot', 'pitchSpeed', 'maxAlt');