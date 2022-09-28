function [] = CAD2pdb(dNode, triad, pdb_path)

% Parameters
n_bp = size(dNode,1);
rot_mat = [0, 0, 1; 1, 0, 0; 0, 1, 0];

% Create strands following the conventions from
% 1. C:\Users\Keyao\Documents\Projects\2012_DNA_Origami\Codes_FileFormats_Alignment\topology2pdb_ver3
% 2. C:\Users\Keyao\Documents\Projects\2012_DNA_Origami\FE_alignment\tools\Local_geometry_v2\transformMat.m
strand(1).isCircular = false;
strand(1).tour = (1 : n_bp)';
strand(1).isMain = true(1, n_bp);
strand(1).seq = repmat('A', 1, n_bp);
for j = 1 : numel(strand(1).tour)
    k = strand(1).tour(j);
    strand(1).R{j} = triad(:,:,k) * rot_mat;
    strand(1).d{j} = dNode(k,:)';
end

strand(2).isCircular = false;
strand(2).tour = (n_bp : -1 : 1)';
strand(2).isMain = false(1, n_bp);
strand(2).seq = repmat('T', 1, n_bp);
for j = 1 : numel(strand(2).tour)
    k = strand(2).tour(j);
    strand(2).R{j} = triad(:,:,k) * rot_mat;    % convert 3DNA conversion to CanDo conversion
    strand(2).d{j} = dNode(k,:)';
end

% Write the PDB file
pdbFinal = pdbGenerate(strand);
pdbFinal = pdbModify_multimodel(pdbFinal);
mypdbwrite_v2(pdb_path, pdbFinal);

end