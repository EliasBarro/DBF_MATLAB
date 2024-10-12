
fileBase = 'C:\Users\austi\Documents\WUDBF_2020\MATLAB';

fileEnd = {...
    %'\wu_proplsn_sim\proplsn_init_scripts\init_proplsn_scorpion_850Kv.m'
    %'\wu_proplsn_sim\proplsn_init_scripts\init_proplsn_sunnysky_380Kv.m'
    %'\wu_proplsn_sim\proplsn_init_scripts\init_proplsn_sunnysky_380Kv_3_blade.m'
    %'\wu_proplsn_sim\proplsn_init_scripts\init_proplsn_biom.m'
    %'\wu_proplsn_sim\proplsn_init_scripts\init_proplsn_hyperion.m'
%     '\wu_proplsn_sim\proplsn_init_scripts\v0\init_proplsn_v0.m'
%     '\wu_proplsn_sim\proplsn_init_scripts\v1\init_proplsn_v1_M2.m'
%     '\wu_proplsn_sim\proplsn_init_scripts\v1\init_proplsn_v1_M3.m'
%     '\wu_proplsn_sim\proplsn_init_scripts\v2\init_proplsn_v2_M2.m'
%     '\wu_proplsn_sim\proplsn_init_scripts\v2\init_proplsn_v2_M3.m'
%     '\wu_proplsn_sim\proplsn_init_scripts\v3\init_proplsn_v3_M3.m'
    '\wu_proplsn_sim\proplsn_init_scripts\v4\init_proplsn_v4_M2.m'
    '\wu_proplsn_sim\proplsn_init_scripts\v4\init_proplsn_v4_M3.m'
    '\wu_proplsn_sim\proplsn_init_scripts\v4\init_proplsn_v4_Flight1a.m'
    '\wu_proplsn_sim\proplsn_init_scripts\v4\init_proplsn_v4_Flight1b.m'
    %'\wu_proplsn_sim\proplsn_init_scripts\init_proplsn_scorpion_330Kv.m'  
    %'\wu_proplsn_sim\proplsn_init_scripts\init_proplsn_scorpion_380Kv.m'
    };

filelist = strcat(fileBase,fileEnd);
% loop through files, generate tables, save off
for indx=1:length(filelist)
    clear proplsn;
    run(filelist{indx});
    try
        prop_perf = get_proplsn_table;
    catch
        disp('Could not find file. Skipping this file.');
        continue;
    end
    [~,basename] = fileparts(filelist{indx});
    savefile = fullfile(wu_km_path,'database',['prop_perf_', basename(14:end)]);
    save(savefile,'prop_perf');
end