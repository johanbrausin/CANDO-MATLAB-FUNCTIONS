clear
clc
addpath C:\Users\User\OneDrive\Escritorio\JOHAN\CANADA-PROJECT\CANDO\CanDo-FRET_v1.0\AV_v2.0
pdbFilename = 'SQUARE-01_DX_52bp_16_cndo_format.pdb';
%#1.1/A:1572@OP1
attAtom.chainSerNo = 1;
attAtom.resSeq = 1572;
attAtom.AtomName = 'OP1';
linkerGeometry(1).L_link = 7;
linkerGeometry(1).w_link = 1;
linkerGeometry(1).R_dye = 3.5;
pdbStruct = pdb2struct_multimodel(pdbFilename);
[posGrid, AV_point] = AV(pdbStruct, attAtom, linkerGeometry);
s = size(AV_point);
x = [];
y = [];
z = [];
for i = 1:s(1)
    x(i) = AV_point(i,1);
    y(i) = AV_point(i,2);
    z(i) = AV_point(i,3);
end
figure;
pcshow([x(:),y(:),z(:)]);
title('Sphere with Default Color Map');
xlabel('X');
ylabel('Y');
zlabel('Z');
%print(posGrid)