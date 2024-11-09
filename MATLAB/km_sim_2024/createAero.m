function aerofinal = createAero(aero_aoa_vs_cl_vs_cd,parasitic_drag)

cl = importdata(aero_aoa_vs_cl_vs_cd);
clhold = cl.data;
climport = clhold(:,[1,16]); % Change columns based on XFLR csv columns (should
                            % be something like RegularFlatWingsStatic)
clfinal = rmmissing(climport);

for i = 1:size(clfinal,1)
    clfinal(i,1) = deg2rad(clfinal(i,1));
end

cd = importdata(aero_aoa_vs_cl_vs_cd);
cdhold = cd.data;
cdimport = cdhold(:,4); % Change columns based on XFLR csv columns (should
                        % be something like RegularFlatWingsStatic)
cdfinal = rmmissing(cdimport);

aerofinal = [clfinal cdfinal];

aerofinal(:,3)=aerofinal(:,3)+parasitic_drag;

end