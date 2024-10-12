<<<<<<< Updated upstream
function aerofinal = createAero(aero_aoa_vs_cl_vs_cd,parasitic_drag)

cl = importdata(aero_aoa_vs_cl_vs_cd);
clhold = cl.data;
climport = clhold(:,[1,3]); % Change columns based on XFLR csv columns (should
                            % be something like RegularFlatWingsStatic)
clfinal = rmmissing(climport);

=======
function aerofinal = createAero(aoa_vs_cl_cd, parasitic_drag)

data = importdata(aoa_vs_cl_cd);
data_hold = data.data;

climport = data_hold(:,[1,7]); % Change columns based on XFLR csv columns (should
                            % be something like BothFlapsandAileronsDown)
clfinal = rmmissing(climport);
>>>>>>> Stashed changes
for i = 1:size(clfinal,1)
    clfinal(i,1) = deg2rad(clfinal(i,1));
end

<<<<<<< Updated upstream
cd = importdata(aero_aoa_vs_cl_vs_cd);
cdhold = cd.data;
cdimport = cdhold(:,4); % Change columns based on XFLR csv columns (should
                        % be something like RegularFlatWingsStatic)
cdfinal = rmmissing(cdimport);

aerofinal = [clfinal cdfinal];

=======
cd = importdata(aoa_vs_cl_cd);
cdhold = cd.data;
cdimport = cdhold(:,8); % Change columns based on XFLR csv columns (should
                        % be something like BothFlapsandAileronsDown)
cdfinal = rmmissing(cdimport);

aerofinal = [clfinal cdfinal];
>>>>>>> Stashed changes
aerofinal(:,3)=aerofinal(:,3)+parasitic_drag;

end