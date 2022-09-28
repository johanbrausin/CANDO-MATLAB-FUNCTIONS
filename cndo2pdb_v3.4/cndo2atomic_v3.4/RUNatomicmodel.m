clear
clc
addpath 'C:\Users\Administrator\Desktop\CANDO-FUNCTIONS\cndo2pdb_v3.4\cndo2atomic_v3.4'; %%%
addpath 'C:\Users\Administrator\Desktop\CANDO-FUNCTIONS\cndo2pdb_v3.4\topology2pdb_ver3_bulge';
param.StrandColor = [190 190 190; 86 180, 233];
%param.StrandColor = [204 121 167; 86 180, 233];   %%OTHERWISAE
param.L_thres = 17415;
param.molmapResolution = 3;
param.WindowSize = [640 480];
param.chimeraEXE = '"C:\ProgramData\Microsoft\Windows\Start Menu\Programs\UCSF ChimeraX 1.4\ChimeraX 1.4.lnk"';
param.chimeraOPTION = '––silent ––script';
CAD_path = 'C:\Users\Administrator\Desktop\CANDO-FUNCTIONS\cndo2pdb_v3.4\cndo2atomic_v3.4\06_Cubeocta_6HB_42bp_Flat_16_CNDO.cndo';
work_DIR = '06_Cubeocta_6HB_42bp_Flat_16_CNDO';
main_cndo2pdb(CAD_path,work_DIR,param)