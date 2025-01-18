function updateSurface(file_name, txt_name)

    loadPlane = file_name; 
    savePlane = file_name;
    filename = txt_name;
    load(loadPlane);
    plane.surfaces = importdata(filename);
    save(savePlane);
   
end