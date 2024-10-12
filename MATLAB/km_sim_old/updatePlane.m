function updatePlane(file_name, txt_name)

    load(file_name, 'plane') %'plane' represent plane struct
    textFile = importdata(txt_name); %open txt file
    plane.name = textFile(1, 1);
    save(file_name, 'plane');

end