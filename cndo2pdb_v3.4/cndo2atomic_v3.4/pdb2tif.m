function [] = pdb2tif(pdb_path, bodyFN, strand, sysParam)

L_thres = sysParam.L_thres;  % staple: L < L_thres

work_dir = fileparts(pdb_path);
tif_1_path = fullfile(work_dir, strcat(bodyFN, '.tif'));
tif_2_path = fullfile(work_dir, strcat(bodyFN, '_x90.tif'));
tif_3_path = fullfile(work_dir, strcat(bodyFN, '_y90.tif'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the UCSF Chimera script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chimeraScr = fullfile(work_dir, strcat(bodyFN, '_chimera.py'));
fid = fopen(chimeraScr, 'w');
% Import the Python interface
fprintf(fid,'from chimera import runCommand\n');

% Open the PDB file
fprintf(fid, 'runCommand(''open %s'')\n', strrep(pdb_path,'\','/'));

% Set the environment
fprintf(fid, 'runCommand(''windowsize %d %d'')\n', sysParam.WindowSize(1), sysParam.WindowSize(2));
fprintf(fid, 'runCommand(''preset apply publication 3'')\n');
fprintf(fid, 'runCommand(''window'')\n');
fprintf(fid, 'runCommand(''scale 0.8'')\n');

% Turn off the original rendering
fprintf(fid, 'runCommand(''~ribbon'')\n');
fprintf(fid, 'runCommand(''~display'')\n');

% Use the new rendering
RGB_scaf = sysParam.StrandColor(1,:)/255;
RGB_stap = sysParam.StrandColor(2,:)/255;
for i = 1:numel(strand)
    if(numel(strand(i).tour) >= L_thres)
        RGB = RGB_scaf;
    else
        RGB = RGB_stap;
    end
    fprintf(fid, 'runCommand(''molmap #0.%d %d'')\n', i, sysParam.molmapResolution);
    fprintf(fid, 'runCommand(''volume #0.%d color %f,%f,%f step 1'')\n', i, RGB(1), RGB(2), RGB(3));
end

% Save as .tif files
fprintf(fid, 'runCommand(''window'')\n');
fprintf(fid, 'runCommand(''scale 0.8'')\n');
fprintf(fid, 'runCommand(''wait'')\n');

fprintf(fid, 'runCommand(''copy file %s tiff dpi 300 supersample 3'')\n', strrep(tif_1_path,'\','/'));
fprintf(fid, 'runCommand(''wait'')\n');

fprintf(fid, 'runCommand(''turn x 90'')\n');
fprintf(fid, 'runCommand(''copy file %s png supersample 3'')\n', strrep(tif_2_path,'\','/'));
fprintf(fid, 'runCommand(''wait'')\n');

fprintf(fid, 'runCommand(''turn x -90'')\n');
fprintf(fid, 'runCommand(''turn y 90'')\n');
fprintf(fid, 'runCommand(''copy file %s png supersample 3'')\n', strrep(tif_3_path,'\','/'));
fprintf(fid, 'runCommand(''wait'')\n');

fprintf(fid, 'runCommand(''turn y -90'')\n');
fprintf(fid, 'runCommand(''close all'')\n');
fprintf(fid, 'runCommand(''stop yes'')\n');
fclose(fid);

runChimera = sprintf('%s %s %s',sysParam.chimeraEXE, sysParam.chimeraOPTION, chimeraScr);
system(runChimera);
    
end
