function [] = render_CAD(filename, dNode, triad, scaf_xover, stap_xover)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constant parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R_node = 0.8;       % [Angstrom] radius of a node
R_e1 = [0.3, 0.6];  % [Angstrom] radii of node axes e1 & e2
R_e3 = [0.4, 1.2];  % [Angstrom] radii of node axis e3
L_e1 = 8.5;         % [Angstrom] lengths of node axes e1 & e2
L_e3 = 2.5;         % [Angstrom] lengths of node axis e3
R_xover = 0.3;      % [Angstrom] radii of crossovers

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Render the FE nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(filename, 'w');

for i = 1 : size(dNode, 1)
    e0 = dNode(i,:)';
    e1 = e0 + triad(:,1,i) * L_e1;
    e2 = e0 + triad(:,2,i) * L_e1;
    e3 = e0 + triad(:,3,i) * L_e3;
    
    % Positions
    fprintf(fid, '.color tan\n');
    fprintf(fid, '.sphere\t%f\t%f\t%f\t%f\n', e0(1), e0(2), e0(3), R_node);
    
    % Orientations
    fprintf(fid, '.color salmon\n');
    fprintf(fid, '.arrow\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n', ...
        e0(1), e0(2), e0(3), e1(1), e1(2), e1(3), R_e1(1), R_e1(2));
    fprintf(fid, '.color sea green\n');
    fprintf(fid, '.arrow\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n', ...
        e0(1), e0(2), e0(3), e2(1), e2(2), e2(3), R_e1(1), R_e1(2));
    fprintf(fid, '.color steel blue\n');
    fprintf(fid, '.arrow\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n', ...
        e0(1), e0(2), e0(3), e3(1), e3(2), e3(3), R_e3(1), R_e3(2));
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Render the crossovers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
assert(size(scaf_xover,2) == 2 && size(stap_xover,2) == 2);

% Scaffold crossovers
fprintf(fid, '.color rosy brown\n');
for i = 1 : size(scaf_xover,1)
    i1 = scaf_xover(i,1);
    i2 = scaf_xover(i,2);
    fprintf(fid, '.cylinder\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n', ...
            dNode(i1,1), dNode(i1,2), dNode(i1,3), dNode(i2,1), dNode(i2,2), dNode(i2,3), R_xover);
end

% Staple crossovers
fprintf(fid, '.color sea green\n');
for i = 1 : size(stap_xover,1)
    i1 = stap_xover(i,1);
    i2 = stap_xover(i,2);
    fprintf(fid, '.cylinder\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n', ...
            dNode(i1,1), dNode(i1,2), dNode(i1,3), dNode(i2,1), dNode(i2,2), dNode(i2,3), R_xover);
end

fclose(fid);

end