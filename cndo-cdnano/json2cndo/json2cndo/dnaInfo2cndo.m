function [] = dnaInfo2cndo(dnaInfo, cndoFN)

fid = fopen(cndoFN, 'w');

% File header
fprintf(fid, '"CanDo (.cndo) file format version 1.0, Keyao Pan, Laboratory for Computational Biology and Biophysics, Massachusetts Institute of Technology, November 2015"\n');
fprintf(fid, '\n');

% dnaInfo.dnaTop
fprintf(fid, 'dnaTop,id,up,down,across,seq\n');
for i = 1 : numel(dnaInfo.dnaTop)
    fprintf(fid, '%d,%d,%d,%d,%d,%s\n', ...
                 i, ...
                 dnaInfo.dnaTop(i).id, ...
                 dnaInfo.dnaTop(i).up, ...
                 dnaInfo.dnaTop(i).down, ...
                 dnaInfo.dnaTop(i).across, ...
                 dnaInfo.dnaTop(i).seq);
end
fprintf(fid, '\n');

% dnaInfo.dnaGeom.dNode
fprintf(fid, 'dNode,"e0(1)","e0(2)","e0(3)"\n');
for i = 1 : size(dnaInfo.dnaGeom.dNode, 1)
    fprintf(fid, '%d,%f,%f,%f\n', ...
                 i, ...
                 dnaInfo.dnaGeom.dNode(i,1), ...
                 dnaInfo.dnaGeom.dNode(i,2), ...
                 dnaInfo.dnaGeom.dNode(i,3));
end
fprintf(fid, '\n');

% dnaInfo.dnaGeom.triad
fprintf(fid, 'triad,"e1(1)","e1(2)","e1(3)","e2(1)","e2(2)","e2(3)","e3(1)","e3(2)","e3(3)"\n');
for i = 1 : size(dnaInfo.dnaGeom.triad, 3)
    fprintf(fid, '%d,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', ...
                 i, ...
                 -dnaInfo.dnaGeom.triad(1,1,i), ...
                 -dnaInfo.dnaGeom.triad(2,1,i), ...
                 -dnaInfo.dnaGeom.triad(3,1,i), ...
                 dnaInfo.dnaGeom.triad(1,2,i), ...
                 dnaInfo.dnaGeom.triad(2,2,i), ...
                 dnaInfo.dnaGeom.triad(3,2,i), ...
                 -dnaInfo.dnaGeom.triad(1,3,i), ...
                 -dnaInfo.dnaGeom.triad(2,3,i), ...
                 -dnaInfo.dnaGeom.triad(3,3,i));
end
fprintf(fid, '\n');

% dnaInfo.dnaGeom.id_nt
fprintf(fid, 'id_nt,id1,id2\n');
for i = 1 : size(dnaInfo.dnaGeom.id_nt, 1)
    fprintf(fid, '%d,%d,%d\n', ...
                 i, ...
                 dnaInfo.dnaGeom.id_nt(i,1), ...
                 dnaInfo.dnaGeom.id_nt(i,2));
end

fclose(fid);

end