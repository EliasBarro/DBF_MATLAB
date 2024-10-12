function [kmpath] = km_path()
%KM_PATH Returns the file path local to the user's computer for
%../sim_dev/km_sim/

fullpath = fileparts(mfilename('fullpath'));
inds = strfind(fullpath,filesep);
kmpathsep = inds(end)-1;            % Go up one from utilities
kmpath = fullpath(1:kmpathsep);
end

