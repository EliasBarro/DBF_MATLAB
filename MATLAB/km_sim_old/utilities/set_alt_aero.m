function set_alt_aero
%SET_ALT_AERO swaps plane.alt_aero table into the primary plane.aero table
%
% intended for flaps, gear down, etc.

% load current plane
plane = current_plane();

% set alternate aero as main aero table
if isfield(plane,'alt_aero')
    plane.aero = plane.alt_aero;
else
    warning('no alternate aero table found in plane structure');
end

% save current plane
current_plane(plane);

end