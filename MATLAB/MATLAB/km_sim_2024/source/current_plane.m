function [plane_out] = current_plane(plane_in)
%CURRENT_PLANE (km_sim toolbox) is mechanism for a pseudo-global variable in 
% the form of a function. current_plane stores a persistent, settable, 
% "current plane" structure that is used for km_sim analysis.
% 
%	SYNTAX: retrive the plane - [plane_out] = current_plane();
%           set the plane     -               current_plane(plane_in);
% 
% see also execute_mission, mission_phase, set_fault

persistent plane_stored;

if nargin==0
    %retriving plane
    plane_out = plane_stored;
elseif nargin==1
    %storing plane
    plane_stored = plane_in;
    plane_out = plane_stored;
else
    error('current_plane called with incorrect inputs')
end

end

