clear
clc
addpath 'C:\Users\Administrator\Desktop\CANDO-FUNCTIONS\cndo2pdb_v3.4\cndo2atomic_v3.4'; %%%%% IN QUOTES ENTER THE DIRECTION OF THE FILE cndo2atomic_v3.4
addpath 'C:\Users\Administrator\Desktop\CANDO-FUNCTIONS\cndo2pdb_v3.4\topology2pdb_ver3_bulge'; %%%%% IN QUOTES ENTER THE DIRECTION OF THE FILE topology2pdb_ver3_bulge

%%%Assign colors for the DNA strands. In this example, if a strand contains less than 100 nucleotides, then a color with the RGB value (190, 190, 190) 
%%%is assigned to this strand. Otherwise a color with the RGB value (204, 121, 167) is assigned.
param.StrandColor = [190 190 190; 86 180, 233]; %%%% 
%param.StrandColor = [204 121 167; 86 180, 233];   %%OTHERWISE
param.L_thres = 100; %%usually works with this number, if don´t, go to the number of nucleotides of your .txt output file of PERDIX or TALOS
param.molmapResolution = 3;
param.WindowSize = [640 480];
param.chimeraEXE = '"C:\ProgramData\Microsoft\Windows\Start Menu\Programs\UCSF ChimeraX 1.4\ChimeraX 1.4.lnk"'; %%%%IN QUOTES ADD THE DIRECTION FILE OF CHIMERAX 
param.chimeraOPTION = '––silent ––script';
%%% ADD THE DIRECTION FILE OF YOUR CNDO FILE AND MAKE SURE THAT YOU´RE UBICATED IN THE "CURRENT FOLDER" WHERE THIS FILE SUPPOSE TO BE 
CAD_path = 'C:\Users\Administrator\Desktop\CANDO-FUNCTIONS\cndo2pdb_v3.4\cndo2atomic_v3.4\06_Cubeocta_6HB_42bp_Flat_16_CNDO.cndo'; 

%%% IN QUOTES ADD A NAME OF THE FOLDER THAT WILL BE CREATED (IF THIS IS MOST THAN YOUR FIRST ATTEMPT RUNING THE CODE, BE SURE TO DELETE THE UNSUCCESFULL FOLDER)
work_DIR = '06_Cubeocta_6HB_42bp_Flat_16_CNDO';
main_cndo2pdb(CAD_path,work_DIR,param)
