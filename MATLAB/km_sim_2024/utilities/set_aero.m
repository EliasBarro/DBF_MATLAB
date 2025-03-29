function set_aero
%SET_AERO swaps the original plane.aero table into the primary plane.aero table
%SET_ALT_AERO must be used beforehand
%
% intended for no flaps, gear up, etc.

% load current plane
plane = current_plane();

% set alternate aero as main aero table
if isfield(plane,'alt_aero')
    plane.aero = plane.old_aero;
else
    warning('no alternate aero table found in plane structure');
end

% save current plane
current_plane(plane);

end