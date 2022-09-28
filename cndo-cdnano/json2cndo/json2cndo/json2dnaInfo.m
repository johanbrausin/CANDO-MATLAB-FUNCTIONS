% Convert a JSON file to a DNA topology data structure
function dnaInfo = json2dnaInfo(jsonPATH, csvPATH, latticeType)
%p_json: http://www.mathworks.com/matlabcentral/fileexchange/25713

% Open and read the JSON file
fid = fopen(jsonPATH);
inString = fscanf(fid,'%c'); 
fclose(fid);
data = p_json(inString);

% Convert the raw data into a master table
% Column 1. UID of the base
% Column 2-4. the current base
% Column 5-7. the 5'-neighbor
% Column 8-10. the 3'-neighbor
% Column 11-13. the Watson-Crick neighbor
% Column 14-16. the Cartesian coordinate
% Column 17. deletion
% Column 18. insertion
[topTable, dNode, triad, id_nt_0] = jsonParser(data,latticeType);
nBase = size(topTable,1);
topTable(:,1) = (1:nBase)';

% Initialize the structure for DNA topology
dnaTop = struct('id',     -1, ...      % uid
                'up',     -1, ...      % uid of the 5' neighbor
                'down',   -1, ...      % uid of the 3' neighbor
                'across', -1, ...      % uid of the Watson-Crick neighbor
                'seq',    'N', ...     % nucleotide identity
                'h',      -1, ...      % helix ID (start from 0)
                'p',      -1, ...      % position in a helix (start from 0)
                'isScaf', false, ...   % Whether this base is on the scaffold
                'skip',   0, ...       % deletion
                'loop',   0);          % insertion
dnaTop = repmat(dnaTop, 1, nBase);

% Create the map table (id_nt)
n_bp = size(dNode,1);
assert(n_bp == size(triad,3) && n_bp == size(id_nt_0,1));
id_nt = zeros(n_bp, 2);
for i = 1 : n_bp
    % Scaffold nucleotide
    tmp = find_row(id_nt_0(i,1:3), topTable(:,2:4));
    assert(numel(tmp)==1);
    id_nt(i,1) = tmp;
    % Staple nucleotide
    tmp = find_row(id_nt_0(i,4:6), topTable(:,2:4));
    assert(numel(tmp)==1);
    id_nt(i,2) = tmp;
end

% Convert the master table into the structure for DNA topology
for i = 1:nBase
    % UID
    dnaTop(i).id = topTable(i,1);
    
    % 5'-neighbor
    tmp = find_row(topTable(i,5:7), topTable(:,2:4));
    if(numel(tmp)==1)
        dnaTop(i).up = tmp;
    elseif(numel(tmp)==0)
        dnaTop(i).up = -1;
    else
        error('Duplication.');
    end
    
    % 3'-neighbor
    tmp = find_row(topTable(i,8:10), topTable(:,2:4));
    if(numel(tmp)==1)
        dnaTop(i).down = tmp;
    elseif(numel(tmp)==0)
        dnaTop(i).down = -1;
    else
        error('Duplication.');
    end
    
    % Watson-Crick neighbor
    tmp = find_row(topTable(i,11:13), topTable(:,2:4));
    if(numel(tmp)==1)
        dnaTop(i).across = tmp;
    elseif(numel(tmp)==0)
        dnaTop(i).across = -1;
    else
        error('Duplication.');
    end
    
    % Helix, position, scaffold/staple
    dnaTop(i).h = topTable(i,2);
    dnaTop(i).p = topTable(i,3);
    if(topTable(i,4) == 0)
        dnaTop(i).isScaf = true;
    elseif(topTable(i,4) == 1)
        dnaTop(i).isScaf = false;
    else
        error('Exception.');
    end
    
    % Deletion
    dnaTop(i).skip = topTable(i,17);
    
    % Insertion
    dnaTop(i).loop = topTable(i,18);
end

% Implement deletions
[dnaTop, dNode, triad, id_nt] = baseDeletion(dnaTop, dNode, triad, id_nt);

% Implement insertions
[dnaTop, dNode, triad, id_nt] = baseInsertion(dnaTop, dNode, triad, id_nt);

% Verification
idList = zeros(size(dnaTop));
for i = 1:numel(dnaTop);
    idList(i) = dnaTop(i).id;
end
assert(numel(unique(idList)) == numel(dnaTop));

% Generate strands
[dnaTop,strand] = buildStrand(dnaTop);
stap_ends = zeros(0,5);
for i = 1 : numel(strand)
    isScaf_strand = dnaTop(strand(i).tour(1)).isScaf;
    for j = 1 : numel(strand(i).tour)
        isScaf_base = dnaTop(strand(i).tour(j)).isScaf;
        assert(~xor(isScaf_strand, isScaf_base));
    end
    
    if(~isScaf_strand)
        assert(~strand(i).isCircular);  % all staples are linear
        h0 = dnaTop(strand(i).tour(1)).h;
        p0 = dnaTop(strand(i).tour(1)).p;
        h1 = dnaTop(strand(i).tour(end)).h;
        p1 = dnaTop(strand(i).tour(end)).p;
        stap_ends = cat(1, stap_ends, [i, h0, p0, h1, p1]);
    end
end

% Assign the sequences
if(~isempty(csvPATH))
    seqInfo = csv2seq(csvPATH);
    for i = 1 : numel(seqInfo)
        tmp = find_row(seqInfo(i).start, stap_ends(:,2:3));
        assert(numel(tmp) == 1);
        assert(norm(seqInfo(i).end - stap_ends(tmp,4:5)) < 1e-10);
        iStrand = stap_ends(tmp,1);
        
        assert(seqInfo(i).len == numel(strand(iStrand).tour));
        for j = 1 : numel(strand(iStrand).tour)
            k = strand(iStrand).tour(j);
            dnaTop(k).seq = seqInfo(i).sequence(j);
            if(dnaTop(k).across >= 0)
                k_across = dnaTop(k).across;
                assert(strcmp(dnaTop(k_across).seq, 'N'));
                dnaTop(k_across).seq = wspair(seqInfo(i).sequence(j));
            end
        end
    end
end

% Output the results
dnaInfo.dnaTop = rmfield(dnaTop, {'h', 'p', 'isScaf', 'skip', 'loop', 'strand', 'residue'});
dnaInfo.dnaGeom.dNode = dNode;
dnaInfo.dnaGeom.triad = triad;
dnaInfo.dnaGeom.id_nt = id_nt;

end


% Parse the JSON data (all helices)
function [topTable, dNode, triad, id_nt_0] = jsonParser(data,latticeType)

% Make sure that there are only two fields: 'name' and 'vstrands'
% assert(numel(fieldnames(data)) == 2);
% assert(isfield(data,'name'));
assert(isfield(data,'vstrands'));

% Read in the topological information in 'vstrands'
topTable = [];
dNode = [];
triad = [];
id_nt_0 = [];
for i = 1:numel(data.vstrands)
    [topTable_0, dNode_0, triad_0, id] = jsonParser_singleHelix(data.vstrands{i},latticeType);
    topTable = cat(1, topTable, topTable_0);
    dNode = cat(1, dNode, dNode_0);
    triad = cat(3, triad, triad_0);
    id_nt_0 = cat(1, id_nt_0, id);
end
assert(numel(topTable) > 0);

end


% Parse the JSON data (single helix)
function [topTable,dNode,triad,id_nt_0] = jsonParser_singleHelix(vstrand,latticeType)

% Check the data format
% assert(numel(fieldnames(vstrand)) == 10);

assert(isfield(vstrand,'num'));
assert(isfield(vstrand,'row'));
assert(isfield(vstrand,'col'));
assert(isfield(vstrand,'scaf'));
assert(isfield(vstrand,'stap'));
assert(isfield(vstrand,'skip'));
assert(isfield(vstrand,'loop'));
% assert(isfield(vstrand,'stap_colors'));
% assert(isfield(vstrand,'scafLoop'));
% assert(isfield(vstrand,'stapLoop'));

% assert(isempty(vstrand.scafLoop));
% assert(isempty(vstrand.stapLoop));

assert(max(cell2mat(vstrand.skip)) <= 0);
assert(min(cell2mat(vstrand.loop)) >= 0);

% Get the maximum number of bases
nLattice = numel(vstrand.scaf);
assert(nLattice == numel(vstrand.stap));
assert(nLattice == numel(vstrand.skip));
assert(nLattice == numel(vstrand.loop));

topTable = zeros(nLattice*2, 18);
dNode = zeros(0,3);
triad = zeros(3,3,0);
id_nt_0 = zeros(0,6);

% Calculate the Cartesian coordinate for each base
[scafXYZ,stapXYZ,dNode_0,triad_0] = initCoor(latticeType, vstrand.row, vstrand.col, vstrand.num, nLattice);
dNode_0 = dNode_0 * 10;     % convert nm to Angstrom (scafXYZ & stapXYZ are still in nm)

% Go through every lattice point
uid = 0;        % UID of the base
for i = 1:nLattice
    currScaf = cell2mat(vstrand.scaf{i});
    currStap = cell2mat(vstrand.stap{i});
    assert((currScaf(1)>=0 && currScaf(2)>=0) || (currScaf(1)<0 && currScaf(2)<0));
    assert((currScaf(3)>=0 && currScaf(4)>=0) || (currScaf(3)<0 && currScaf(4)<0));
    assert((currStap(1)>=0 && currStap(2)>=0) || (currStap(1)<0 && currStap(2)<0));
    assert((currStap(3)>=0 && currStap(4)>=0) || (currStap(3)<0 && currStap(4)<0));
    
    if(currScaf(1)>=0 || currScaf(3)>=0)
        % The base exists in the scaffold strand
        uid = uid+1;
        topTable(uid,1) = uid;
        topTable(uid,2:4) = [vstrand.num, i-1, 0];    % helix ID, lattice ID, 0 for scaffold, 1 for staple
        uid_currScaf = topTable(uid,2:4);
        
        % The 5'-neighbor
        topTable(uid,5:7) = [currScaf(1:2), 0];
        
        % The 3'-neighbor
        topTable(uid,8:10) = [currScaf(3:4), 0];
        
        % The Watson-Crick neighbor
        if(currStap(1)>=0 || currStap(3)>=0)
            topTable(uid,11:13) = [vstrand.num, i-1, 1];
        else
            topTable(uid,11:13) = [-1, -1, 1];
        end
        
        % The Cartesian coordinate
        topTable(uid,14:16) = scafXYZ(i,:);
        
        % The deletion
        topTable(uid,17) = vstrand.skip{i};
        assert(topTable(uid,17)==0 || topTable(uid,17)==-1);
        
        % The insertion
        topTable(uid,18) = vstrand.loop{i};
        assert(topTable(uid,18)>=0);
    end
    
    if(currStap(1)>=0 || currStap(3)>=0)
        % The base exists in the staple strand
        uid = uid+1;
        topTable(uid,1) = uid;
        topTable(uid,2:4) = [vstrand.num, i-1, 1];
        uid_currStap = topTable(uid,2:4);
        
        % The 5'-neighbor
        topTable(uid,5:7) = [currStap(1:2), 1];
        
        % The 3'-neighbor
        topTable(uid,8:10) = [currStap(3:4), 1];
        
        % The Watson-Crick neighbor
        if(currScaf(1)>=0 || currScaf(3)>=0)
            topTable(uid,11:13) = [vstrand.num, i-1, 0];
        else
            topTable(uid,11:13) = [-1, -1, 0];
        end
        
        % The Cartesian coordinate
        topTable(uid,14:16) = stapXYZ(i,:);
        
        % The deletion
        topTable(uid,17) = vstrand.skip{i};
        assert(topTable(uid,17)==0 || topTable(uid,17)==-1);
        
        % The insertion
        topTable(uid,18) = vstrand.loop{i};
        assert(topTable(uid,18)>=0);
    end
    
    % Check if a basepair exists
    if((currScaf(1)>=0 || currScaf(3)>=0) && (currStap(1)>=0 || currStap(3)>=0))
        dNode = cat(1, dNode, dNode_0(i,:));
        triad = cat(3, triad, triad_0(:,:,i));
        id_nt_0 = cat(1, id_nt_0, [uid_currScaf, uid_currStap]);
    end
end

topTable(uid+1:end,:) = [];

end


% Find a row vector in all the rows in a matrix
function ind = find_row(a,M)

assert(size(a,1)==1 && size(a,2)==size(M,2));
tmp = bsxfun(@minus,M,a);
s = sum(abs(tmp),2);
ind = find(s==0);

end


% Calculate the Cartesian coordinate of a base
function [scafXYZ,stapXYZ,dNode,triad] = initCoor(latticeType, row, col, strandNum, nLattice)

% Five constant parameters
rStrand = 1.25;     % Half the distance between the axes of two neighboring DNA helices
rHelix = 1;         % Radius of DNA helices (nm)
distBP = 0.34;      % Rise between two neighboring base-pairs (nm)
angBP = 360/10.5;   % Twisting angle between two neighboring base-pairs (degree)
angMinor = 120;     % Angle of the minor groove (degree) (http://en.wikibooks.org/wiki/Structural_Biochemistry/Nucleic_Acid/DNA/DNA_structure)

% Positions of the scaffold nucleotide and staple nucleotide
% in the local reference frame
scafLocal = rHelix * [cos(180-angMinor/2), sin(180-angMinor/2), 0]';
stapLocal = rHelix * [cos(180+angMinor/2), sin(180+angMinor/2), 0]';

% See CanDo ver 2.3, adinaTranslate.m
% starting x, y, z coordinates
if strcmp(latticeType,'honeycomb')
    xpos = col*rStrand*sqrt(3);
    zpos = -rStrand*3*row;
    if(((mod(row,2)==0)&&(mod(col,2)==0))||((mod(row,2)~=0)&&(mod(col,2)~=0)))
        zpos = zpos+rStrand;
    end
    if mod(strandNum,2)==0 % even numbered strand
        initAng = -30 + angBP/2;
    else % odd numbered strand
        initAng = 150 + angBP/2;
    end
elseif strcmp(latticeType,'square')
    xpos = col*rStrand*2;
    zpos = -rStrand*2*row;
    if mod(strandNum,2)==0 % even numbered strand
        initAng = 180 + angBP/2;
    else % odd numbered strand
        initAng = 0 + angBP/2;
    end
else
    error('Unknown lattice type!');
end
initCoordAngStrand = [xpos 0.0 zpos initAng];

% XYZ coordinate and twisting angle for each base-pair
assert(nLattice>=0);
scafXYZ = zeros(nLattice,3);
stapXYZ = zeros(nLattice,3);
dNode = zeros(nLattice,3);
triad = zeros(3,3,nLattice);
% angCorrection = (180-angMinor)/2;

if(mod(strandNum,2)==0) % even numbered strand
    e3 = [0, 1, 0]';
else % odd numbered strand
    e3 = [0, -1, 0]';
end

for i = 0:nLattice-1
    refCoord = initCoordAngStrand + [0, distBP*i, 0, angBP*i];
    % Basepair position
    dNode(i+1,:) = refCoord(1:3);
    % Basepair orientation
    e2 = [cos(-deg2rad(refCoord(4))), 0, sin(-deg2rad(refCoord(4)))]';
    e1 = cross(e2, e3);
    triad(:,:,i+1) = [e1, e2, e3];
    % Scaffold & staple nucleotide positions
    scafXYZ(i+1,:) = dNode(i+1,:) + (triad(:,:,i+1)*scafLocal)';
    stapXYZ(i+1,:) = dNode(i+1,:) + (triad(:,:,i+1)*stapLocal)';
    
%     centerXYZ = refCoord(1:3);
%     if(mod(strandNum,2)==0)
%         angScaf = refCoord(4) + angCorrection;
%         angStap = angScaf + angMinor;
%     else
%         angScaf = refCoord(4) - angCorrection;
%         angStap = angScaf - angMinor;
%     end
%     scafXYZ(i+1,:) = centerXYZ + rHelix*[cos(-angScaf/180*pi), 0, sin(-angScaf/180*pi)];
%     stapXYZ(i+1,:) = centerXYZ + rHelix*[cos(-angStap/180*pi), 0, sin(-angStap/180*pi)];
end

end


% Implement deletions
function [dnaTop, dNode, triad, id_nt] = baseDeletion(dnaTop, dNode, triad, id_nt)

nBase = numel(dnaTop);
skip = zeros(nBase,1);
for i = 1:nBase
    skip(i) = dnaTop(i).skip;
end

skipList = [];
skipList_bp = [];
currSkip = find(skip~=0,1);
while(~isempty(currSkip))
    assert(skip(currSkip) == -1);
    
    % Find the four neighboring bases
    neighborUp = dnaTop(currSkip).up;
    neighborDown = dnaTop(currSkip).down;
    currSkipAcross = dnaTop(currSkip).across;
    if(currSkipAcross >= 0)
        neighborAcrossUp = dnaTop(currSkipAcross).up;
        neighborAcrossDown = dnaTop(currSkipAcross).down;
    else
        neighborAcrossUp = -1;
        neighborAcrossDown = -1;
    end
    
    % Change the connectivity
    if(neighborUp >= 0)
        dnaTop(neighborUp).down = neighborDown;
    end
    if(neighborDown >= 0)
        dnaTop(neighborDown).up = neighborUp;
    end
    if(neighborAcrossUp >= 0)
        dnaTop(neighborAcrossUp).down = neighborAcrossDown;
    end
    if(neighborAcrossDown >= 0)
        dnaTop(neighborAcrossDown).up = neighborAcrossUp;
    end
    
    % Add to a list of the bases to be deleted
    if(currSkipAcross >= 0)
        % Update the skip list for bases
        skipList = cat(1, skipList, [currSkip; currSkipAcross]);
        % Update the skip list for basepairs
        [ti, tj] = find(id_nt == currSkip);
        assert(numel(ti)==1 && numel(tj)==1);
        assert(id_nt(ti,3-tj) == currSkipAcross);
        skipList_bp = cat(1, skipList_bp, ti);
        % Remove the currect base from the remaining skips
        skip([currSkip, currSkipAcross]) = 0;
    else
        % Update the skip list for bases
        skipList = cat(1, skipList, currSkip);
        % Remove the currect base from the remaining skips
        skip(currSkip) = 0;
    end
    
    % Find the next base to be deleted
    currSkip = find(skip~=0,1);
end

% Delete bases
dnaTop(skipList) = [];
dNode(skipList_bp,:) = [];
triad(:,:,skipList_bp) = [];
id_nt(skipList_bp,:) = [];

% Change the base IDs
[dnaTop, id_nt] = baseID_squeeze(dnaTop, id_nt);

end


% Implement insertions
function [dnaTop, dNode, triad, id_nt] = baseInsertion(dnaTop, dNode, triad, id_nt)

distBP = 3.4;       % Rise between two neighboring base-pairs (Angstrom)
angBP = 360/10.5;   % Twisting angle between two neighboring base-pairs (degree)

% dy = [0, 0.1, 0]';
nBase = numel(dnaTop);
loop = zeros(nBase,1);
for i = 1:nBase
    loop(i) = dnaTop(i).loop;
end

currLoop = find(loop~=0,1);
while(~isempty(currLoop))
    assert(loop(currLoop) > 0);
    currLoopAcross = dnaTop(currLoop).across;
    if(currLoopAcross >= 0)     % dsDNA
        assert(loop(currLoopAcross) > 0);
        [ti,tj] = find(id_nt == currLoop);
        assert(numel(ti)==1 && numel(tj)==1);
        assert(id_nt(ti,3-tj) == currLoopAcross);
        % Make sure that the 3'-neighbor of the base 'currLoop' is to the
        % right
        tmp_a = triad(:,3,ti);
        if(norm(tmp_a - [0;1;0]) < 1e-10)  % reference axis e3 pointing to the right
            if(tj == 2) % base 'currLoop' is a non-preferred base
                tmp_nt = currLoop;
                currLoop = currLoopAcross;
                currLoopAcross = tmp_nt;
            end
        elseif(norm(tmp_a - [0;-1;0]) < 1e-10) % reference axis e3 pointing to the left
            if(tj == 1) % base 'currLoop' is a preferred base
                tmp_nt = currLoop;
                currLoop = currLoopAcross;
                currLoopAcross = tmp_nt;
            end
        else
            error('Exception.');
        end
        
    end
    
    % Find the four neighboring bases
    neighborDown = dnaTop(currLoop).down;
    
    if(currLoopAcross >= 0)
        neighborAcrossUp = dnaTop(currLoopAcross).up;
    else
        neighborAcrossUp = -1;
    end
    
    % Insert bases
    if(currLoopAcross >= 0)
        nInsert = loop(currLoop)*2;
    else
        nInsert = loop(currLoop);
    end
    
    dnaTop_insert = struct('id',     -1, ...      % uid
                           'up',     -1, ...      % uid of the 5' neighbor
                           'down',   -1, ...      % uid of the 3' neighbor
                           'across', -1, ...      % uid of the Watson-Crick neighbor
                           'seq',    'N', ...     % nucleotide identity
                           'h',      -1, ...      % helix ID (start from 0)
                           'p',      -1, ...      % position in a helix (start from 0)
                           'isScaf', false, ...   % Whether this base is on the scaffold
                           'skip',   0, ...       % deletion
                           'loop',   0);          % insertion
    dnaTop_insert = repmat(dnaTop_insert, 1, nInsert);
    
    % Case I. Insert dsDNA
    if(currLoopAcross >= 0)
        assert(mod(nInsert,2)==0);
        
        % Create basepair information
        [ti,tj] = find(id_nt == currLoop);
        assert(numel(ti)==1 && numel(tj)==1);
        assert(id_nt(ti,3-tj) == currLoopAcross);
        
        % Positions and orientations
        dNode_1 = dNode(ti,:);
        triad_1 = triad(:,:,ti);
        dNode_2 = dNode_1 + [0, distBP, 0];
        triad_2 = vrrotvec2mat([0, 1, 0, deg2rad(angBP)]) * triad_1;
        [dNode_insert, triad_insert] = bp_interp(dNode_1, triad_1, dNode_2, triad_2, nInsert/2);
        % Map table
        id_nt_insert = reshape(nBase+(1:nInsert), 2, nInsert/2)';
        if(norm(triad_1(:,3) - [0;1;0]) < 1e-10)
            % Do nothing
        elseif(norm(triad_1(:,3) - [0;-1;0]) < 1e-10)  % base 'currLoop' is non-preferred
            id_nt_insert = fliplr(id_nt_insert);
        else
            error('Exception.');
        end
        
        % Insert new basepairs
        dNode = cat(1, dNode, dNode_insert);
        triad = cat(3, triad, triad_insert);
        id_nt = cat(1, id_nt, id_nt_insert);

        for i = 1:2:nInsert
            dnaTop_insert(i).id = nBase+i;
            dnaTop_insert(i+1).id = nBase+i+1;
            
            if(i==1)
                dnaTop(currLoop).down = nBase+1;
                dnaTop(currLoopAcross).up = nBase+2;
                dnaTop_insert(i).up = currLoop;
                dnaTop_insert(i+1).down = currLoopAcross;
            else
                assert(i-2 > 0);
                dnaTop_insert(i).up = nBase+i-2;
                dnaTop_insert(i+1).down = nBase+i-1;
            end
            
            if(i==nInsert-1)
                if(neighborDown>=0)
                    dnaTop(neighborDown).up = nBase+nInsert-1;
                end
%                 disp([currLoop,currLoopAcross,neighborAcrossUp])
                if(neighborAcrossUp>=0)
                    dnaTop(neighborAcrossUp).down = nBase+nInsert;
                end
                dnaTop_insert(i).down = neighborDown;
                dnaTop_insert(i+1).up = neighborAcrossUp;
            else
                assert(i+3 <= nInsert);
                dnaTop_insert(i).down = nBase+i+2;
                dnaTop_insert(i+1).up = nBase+i+3;
            end
                
            dnaTop_insert(i).across = dnaTop_insert(i+1).id;
            dnaTop_insert(i+1).across = dnaTop_insert(i).id;
            
            dnaTop_insert(i).h = dnaTop(currLoop).h;
            dnaTop_insert(i+1).h = dnaTop(currLoopAcross).h;
            dnaTop_insert(i).p = dnaTop(currLoop).p;
            dnaTop_insert(i+1).p = dnaTop(currLoopAcross).p;
            dnaTop_insert(i).isScaf = dnaTop(currLoop).isScaf;
            dnaTop_insert(i+1).isScaf = dnaTop(currLoopAcross).isScaf;
            
            dnaTop_insert(i).skip = 0;
            dnaTop_insert(i+1).skip = 0;
            dnaTop_insert(i).loop = 0;
            dnaTop_insert(i+1).loop = 0;
        end
    
    % Case II. Insert ssDNA
    else
        for i = 1:nInsert
            dnaTop_insert(i).id = nBase+i;
            
            if(i==1)
                dnaTop(currLoop).down = nBase+1;
                dnaTop_insert(i).up = currLoop;
            else
                dnaTop_insert(i).up = nBase+i-1;
            end
            
            if(i==nInsert)
                dnaTop(neighborDown).up = nBase+nInsert;
                dnaTop_insert(i).down = neighborDown;
            else
                dnaTop_insert(i).down = nBase+i+1;
            end
            
            dnaTop_insert(i).across = -1;
            
%             pr = i/(nInsert+1);
%             dnaTop_insert(i).xyz = dnaTop(currLoop).xyz*pr + dnaTop(neighborDown).xyz*(1-pr);
%             dnaTop_insert(i).xyz = dnaTop(currLoop).xyz + dy*i;
            
            dnaTop_insert(i).skip = 0;
            dnaTop_insert(i).loop = 0;
        end
    end
    
    % Remove records
    if(currLoopAcross >= 0)
        loop([currLoop, currLoopAcross]) = 0;
    else
        loop(currLoop) = 0;
    end
    
    % Expand the graph
    dnaTop = cat(2, dnaTop, dnaTop_insert);
    nBase = numel(dnaTop);
    
    % Find the next base to be inserted
    currLoop = find(loop~=0,1);
end

end


% Let base ID be 1:numel(dnaTop)
function [dnaTop_new, id_nt_new] = baseID_squeeze(dnaTop, id_nt)

% Convert topology structure into an array
nBase = numel(dnaTop);
connectTable = zeros(nBase,4);
for i = 1:nBase
    connectTable(i,:) = [dnaTop(i).id, dnaTop(i).up, dnaTop(i).down, dnaTop(i).across];
end

% Sort this array
[~,ind] = sort(connectTable(:,1));
connectTable = connectTable(ind,:);

% Modify base IDs
id_nt_new = id_nt;
connectTable_new = -2*ones(size(connectTable));
connectTable_new(connectTable == -1) = -1;
for i = 1:size(connectTable,1)
    connectTable_new(connectTable == connectTable(i,1)) = i;
    id_nt_new(id_nt == connectTable(i,1)) = i;
end
assert(isempty(find(connectTable_new == -2, 1)));

% Save as a new topology structure
dnaTop_new = struct('id',     -1, ...      % uid
                    'up',     -1, ...      % uid of the 5' neighbor
                    'down',   -1, ...      % uid of the 3' neighbor
                    'across', -1, ...      % uid of the Watson-Crick neighbor
                    'seq',    'N', ...     % nucleotide identity
                    'h',      -1, ...      % helix ID (start from 0)
                    'p',      -1, ...      % position in a helix (start from 0)
                    'isScaf', false, ...   % Whether this base is on the scaffold
                    'skip',   0, ...       % deletion
                    'loop',   0);          % insertion
dnaTop_new = repmat(dnaTop_new, 1, nBase);

for i = 1:nBase
    dnaTop_new(i).id = connectTable_new(i,1);
    dnaTop_new(i).up = connectTable_new(i,2);
    dnaTop_new(i).down = connectTable_new(i,3);
    dnaTop_new(i).across = connectTable_new(i,4);
    dnaTop_new(i).h = dnaTop(ind(i)).h;
    dnaTop_new(i).p = dnaTop(ind(i)).p;
    dnaTop_new(i).isScaf = dnaTop(ind(i)).isScaf;
    dnaTop_new(i).skip = dnaTop(ind(i)).skip;
    dnaTop_new(i).loop = dnaTop(ind(i)).loop;
end

end


% Interpolate the positions and orientations
% See the function fit_R_d in C:\Users\Keyao Pan\Documents\Projects\
%     2012_DNA_Origami\Codes_FileFormats_Alignment\topology2pdb_ver3_bulge\
%     generateBulgeDOF.m
function [dNode_interp, triad_interp] = bp_interp(dNode_1, triad_1, dNode_2, triad_2, N)

dNode_interp = zeros(N,3);
triad_interp = zeros(3,3,N);

% Rotation matrix R
% R * triad_1 = triad_2
R = triad_2 / triad_1;
v = vrrotmat2vec(R);
a = v(1:3);
theta = v(4);

% Calculate for dNode_interp and triad_interp
for i = 1 : N
    dNode_interp(i,:) = (dNode_1*(N+1-i) + dNode_2*i) / (N+1);
    triad_interp(:,:,i) = vrrotvec2mat([a, theta*i/(N+1)]) * triad_1;
end

end