%Creates a Plane Struct and stores it into a .mat file as specified by store_file_name

% store_file_name -      string specifying the name of the plane struct
% aero_aoa_vs_cl_vs_cd - csv of angle of attack vs coefficient of lift and 
%                        coefficient of drag WITHOUT PARASITIC DRAG ADDED (from XFLR)
% single_var_csv_file -  CSV file with static variables that are put into plane struct
% thrust_vs_velo -       CSV file containing thrust values in the first column and 
%                        corresponding velocity values in the second column
function createPlaneStruct(store_file_name, aero_aoa_vs_cl_vs_cd, single_var_csv_file, thrust_vs_velo)
    %initialize struct
    plane = struct; 
    
    %initialize single variables
    single_vars = readmatrix(single_var_csv_file);
    plane.Amps = single_vars(1); %amps
    plane.maxAlt = single_vars(2); %meters
    plane.g_limit = single_vars(3); %m/s^2
    plane.empty_W = 9.80665*single_vars(4); %kg
    plane.spec_fuel_W = single_vars(5); %kg
    plane.e_cap = single_vars(6); %milliAmp hours
    plane.weight = 9.80665*single_vars(7); %kg
    plane.S_ref = single_vars(8); %m^2
    plane.parasitic_drag = single_vars(9); %unitless
    
    %initialize aero and alt_aero
    plane.aero = createAero(aero_aoa_vs_cl_vs_cd, plane.parasitic_drag);
    plane.alt_aero = createAltAero(aero_aoa_vs_cl_vs_cd, plane.parasitic_drag);
    
    
    %initialize prop perf
    plane.prop_perf = gen_linear_prop_table(plane.Amps*-1000/3600, plane.maxAlt, thrust_vs_velo);
    
    %initialize static variables
    plane.TO_frict_coeff = 0.0030;

    save(store_file_name, 'plane');
end
