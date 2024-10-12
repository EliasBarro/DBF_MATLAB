function aerofinal = createAero(aoa_vs_cl, aoa_vs_cd,parasitic_drag)

cl = importdata(aoa_vs_cl);
clhold = cl.data;
climport = clhold(:,[1,2]); % Change columns based on XFLR csv columns (should
                            % be something like RegularFlatWingsStatic)
clfinal = rmmissing(climport);

for i = 1:size(clfinal,1)
    clfinal(i,1) = deg2rad(clfinal(i,1));
end

cd = importdata(aoa_vs_cd);
cdhold = cd.data;
cdimport = cdhold(:,2); % Change columns based on XFLR csv columns (should
                        % be something like RegularFlatWingsStatic)
cdfinal = rmmissing(cdimport);

aerofinal = [clfinal cdfinal];

aerofinal(:,3)=aerofinal(:,3)+parasitic_drag;

end