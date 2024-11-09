function aerofinal = createAltAero(aoa_vs_cl, aoa_vs_cd,parasitic_drag)

cl = importdata(aoa_vs_cl);
clhold = cl.data;
climport = clhold(:,[7,8]); % Change columns based on XFLR csv columns (should
                            % be something like BothFlapsandAileronsDown)
clfinal = rmmissing(climport);

for i = 1:size(clfinal,1)
    clfinal(i,1) = deg2rad(clfinal(i,1));
end

cd = importdata(aoa_vs_cd);
cdhold = cd.data;
cdimport = cdhold(:,8); % Change columns based on XFLR csv columns (should
                        % be something like BothFlapsandAileronsDown)
cdfinal = rmmissing(cdimport);

aerofinal = [clfinal cdfinal];

aerofinal(:,3)=aerofinal(:,3)+parasitic_drag;

end