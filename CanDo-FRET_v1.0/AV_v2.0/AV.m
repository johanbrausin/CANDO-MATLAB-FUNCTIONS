%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main function
% Ref. Kalinin 2012, Nature Methods 9, 1218 (See P.10 in SI for details)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [posGrid,pos] = AV(pdbStruct, attachPos, dyeGeom)

% Constant parameters (see Kalinin 2012, SI P.10)
param.AVGrid = 0.2;       % Default spacing (ratio of min(L_link,w_link,R_dye))
param.minGrid = 0.4;      % The lower bound of the grid spacing size (unit: Angstrom)
%param.allowSphere = 0.5;  % See the text
param.allowSphere = 1;
param.searchNode = 3;     % See the text
param.vdwC = 1.70;        % vdw radius (unit: Angstrom)
param.vdwH = 1.09;
param.vdwO = 1.52;
param.vdwN = 1.55;
param.vdwP = 1.80;
param.vdwS = 1.80;

tic;
% Step 1. Build a 3D grid
[posGrid, posGridOrigin, gridSpace, mmXYZR] = build3DGrid(pdbStruct, attachPos, dyeGeom, param);
toc;
% Step 2. Clash check
posGrid = clashCheck(posGrid, posGridOrigin, gridSpace, mmXYZR, dyeGeom, param);
toc;
% Step 3. Find the shortest route distance for each point in the grid
posGrid = routeDist(posGrid, param);
toc;
% Step 4. Find all the positions in the AV
pos = findPos(posGrid, dyeGeom);
toc;

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build a 3D grid (posGrid) with the origin at the attachment atom, a list 
% of macromolecular atom positions/vdw radii (mmXYZR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [posGrid, posGridOrigin, gridSpace, mmXYZR] = build3DGrid(pdbStruct, attachPos, dyeGeom, param)

% Find the origin of the 3D grid
tmp = find(pdbStruct.chainSerNo == attachPos.chainSerNo & ...
           pdbStruct.resSeq == attachPos.resSeq & ...
           strcmp(pdbStruct.AtomName, attachPos.AtomName));
if(numel(tmp) ~= 1)
    error('Non-unique attachment.');
end
posGridOrigin = pdbStruct.XYZ(tmp,:);
       
% Find which macromolecular atoms may clash with the dye
nAtom = size(pdbStruct.XYZ,1);
mmXYZR = zeros(nAtom,5);
for i = 1:nAtom
    currXYZ = pdbStruct.XYZ(i,:);               % XYZ
    currAtomName = pdbStruct.AtomName{i};       % Atom name
    if(strcmp(currAtomName(1), 'C'))            % vdw radius
        currR = param.vdwC;
    elseif(strcmp(currAtomName(1), 'H'))
        currR = param.vdwH;
    elseif(strcmp(currAtomName(1), 'O'))
        currR = param.vdwO;
    elseif(strcmp(currAtomName(1), 'N'))
        currR = param.vdwN;
    elseif(strcmp(currAtomName(1), 'P'))
        currR = param.vdwP;
    elseif(strcmp(currAtomName(1), 'S'))
        currR = param.vdwS;
    elseif((strcmp(currAtomName(1), '1') || strcmp(currAtomName(1), '2')) && strcmp(currAtomName(2), 'H'))
        currR = param.vdwH;
    else
        error('Illegal atom name.');
    end
    
    % Whether the current atom may clash with the dye?
    if(norm(currXYZ - posGridOrigin) <= dyeGeom.L_link + max(dyeGeom.R_dye,dyeGeom.L_link/2) + currR)
        mmXYZR(i,1:3) = currXYZ;
        mmXYZR(i,4) = currR;
        mmXYZR(i,5) = 1;
    end
end

mmXYZR(mmXYZR(:,5)==0, :) = [];
mmXYZR(:,5) = [];

% Build the 3D grid (2N+1)*(2N+1)*(2N+1)
gridSpace = min([dyeGeom.L_link, dyeGeom.w_link, dyeGeom.R_dye] * param.AVGrid);
gridSpace = max(gridSpace, param.minGrid);
N = ceil(dyeGeom.L_link / gridSpace);

posGrid(2*N+1,2*N+1,2*N+1).XYZ = zeros(1,3);
posGrid(2*N+1,2*N+1,2*N+1).isClashLink = 0;
posGrid(2*N+1,2*N+1,2*N+1).isClashDye = 0;
posGrid(2*N+1,2*N+1,2*N+1).isVisit = 0;
posGrid(2*N+1,2*N+1,2*N+1).d_route = Inf;
for i = 1:(2*N+1)
    for j = 1:(2*N+1)
        for k = 1:(2*N+1)
            posGrid(i,j,k).XYZ = posGridOrigin + gridSpace * ([i j k] - N - 1);
            posGrid(i,j,k).isClashLink = 0;
            posGrid(i,j,k).isClashDye = 0;
            posGrid(i,j,k).isVisit = 0;
            posGrid(i,j,k).d_route = Inf;
        end
    end
end

% Set the start point (origin) for the routing algorithm
posGrid(N+1,N+1,N+1).d_route = 0;

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clash check
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function posGrid = clashCheck(posGrid, posGridOrigin, gridSpace, mmXYZR, dyeGeom, param)

% Set the origin of the grid
originIJK = (size(posGrid)+1)/2;

% Set the limits
iMin = 1;
jMin = 1;
kMin = 1;
[iMax,jMax,kMax] = size(posGrid);

% Go through each macromolecular atom
for s = 1:size(mmXYZR,1);
    r0 = mmXYZR(s,4) + max(dyeGeom.R_dye,dyeGeom.L_link);
    i1 = originIJK(1) - ceil(r0/gridSpace);
    i2 = originIJK(1) + ceil(r0/gridSpace);
    j1 = originIJK(2) - ceil(r0/gridSpace);
    j2 = originIJK(2) + ceil(r0/gridSpace);
    k1 = originIJK(3) - ceil(r0/gridSpace);
    k2 = originIJK(3) + ceil(r0/gridSpace);
    
    for i = max(iMin,i1) : min(iMax,i2)
        for j = max(jMin,j1) : min(jMax,j2)
            for k = max(kMin,k1) : min(kMax,k2)
                if(0 == posGrid(i,j,k).isClashLink || 0 == posGrid(i,j,k).isClashDye)
                    r = norm(posGrid(i,j,k).XYZ - mmXYZR(s,1:3));
                    % If macromolecular atom clashs with the linker
                    if(r < dyeGeom.w_link/2 + mmXYZR(s,4) && norm(posGrid(i,j,k).XYZ - posGridOrigin) > param.allowSphere * dyeGeom.w_link)
                        posGrid(i,j,k).isClashLink = 1;
                    end
                    % If macromolecular atom clashs with the dye
                    if(r < dyeGeom.R_dye + mmXYZR(s,4))
                        posGrid(i,j,k).isClashDye = 1;
                    end
                end
                if(posGrid(i,j,k).isClashLink==1 && posGrid(i,j,k).isClashDye == 0)
                    error('debug');
                end
            end
        end
    end
end

%for i = 1:size(posGrid,1)
%    for j = 1:size(posGrid,2)
%        for k = 1:size(posGrid,3)
%            for s = 1:size(mmXYZR,1)
%                r = norm(posGrid(i,j,k).XYZ - mmXYZR(s,1:3));
%                if(r < dyeGeom.R_dye + mmXYZR(s,4))
%                    posGrid(i,j,k).isClash = 1;
%                end
%            end
%            if(norm(posGrid(i,j,k).XYZ - posGridOrigin) < param.allowSphere * dyeGeom.w_link)
%                posGrid(i,j,k).isClash = 0;
%                %posGrid(i,j,k).isVisit = 0;
%                %posGrid(i,j,k).d_route = norm(posGrid(i,j,k).XYZ - posGridOrigin);
%            end
%        end
%    end
%end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the shortest route distance from the attachment point to the center
% of the dye
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function posGrid = routeDist(posGrid, param)

% Set the limits
iMin = 1;
jMin = 1;
kMin = 1;
[iMax,jMax,kMax] = size(posGrid);

% Define the neighbors
currIJK = (param.searchNode+1)*ones(1,3);
neighbor = findNeighbor(posGrid, currIJK, param.searchNode);
neighbor = bsxfun(@minus, neighbor, currIJK);

% Initialization
originIJK = (size(posGrid)+1)/2;
posGrid(originIJK(1),originIJK(2),originIJK(3)).d_route = 0;
queue = zeros(0,3);
queue = enqueue(queue,originIJK);
posGrid(originIJK(1),originIJK(2),originIJK(3)).isVisit = 1;

while(~isempty(queue))
    [queue,currIJK] = dequeue(queue);
    i0 = currIJK(1);
    j0 = currIJK(2);
    k0 = currIJK(3);
    
    for s = 1:size(neighbor,1)
        % [i j k] is a grid point in the neighborhood of [i0 j0 k0].
        i = i0 + neighbor(s,1);
        j = j0 + neighbor(s,2);
        k = k0 + neighbor(s,3);
        if(i>=iMin && i<=iMax && ...
           j>=jMin && j<=jMax && ...
           k>=kMin && k<=kMax && ...
           0 == posGrid(i,j,k).isClashLink) %&& 0 == posGrid(i,j,k).isVisit)
            % If the current neighbor does not clash with the macromolecular atoms and
            % is new (has never entered the queue) ...
            if(0 == posGrid(i,j,k).isVisit)
                queue = enqueue(queue, [i j k]);
                posGrid(i,j,k).isVisit = 1;
            end
            
            d = norm(posGrid(i,j,k).XYZ - posGrid(i0,j0,k0).XYZ);
            %if(posGrid(i,j,k).d_route > posGrid(i0,j0,k0).d_route + d)
            %    posGrid(i,j,k).d_route = posGrid(i0,j0,k0).d_route + d;
            %end
            if(posGrid(i0,j0,k0).d_route > posGrid(i,j,k).d_route + d)
                posGrid(i0,j0,k0).d_route = posGrid(i,j,k).d_route + d;
            end
        end
    end
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the neighbors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pos = findPos(posGrid, dyeGeom)

pos = [];

for i = 1:size(posGrid,1)
    for j  = 1:size(posGrid,2)
        for k = 1:size(posGrid,3)
            if(0 == posGrid(i,j,k).isClashDye && posGrid(i,j,k).d_route < dyeGeom.L_link)
                pos = cat(1, pos, posGrid(i,j,k).XYZ);
            end
        end
    end
end


end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the neighbors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function neighbor = findNeighbor(posGrid, currIJK, depth)

% Constants
move = [1 0 0; ...
        -1 0 0; ...
        0 1 0; ...
        0 -1 0; ...
        0 0 1; ...
        0 0 -1];

% Set the limits
iMin = 1;
jMin = 1;
kMin = 1;
[iMax,jMax,kMax] = size(posGrid);

% Initialize the queue and isVisit
isVisit = false(iMax,jMax,kMax);
queue = currIJK;
isVisit(currIJK(1),currIJK(2),currIJK(3)) = true;
neighbor = zeros(0,3);

% Search for neighbors
for loop = 1:depth
    queueNext = zeros(0,3);
    while(~isempty(queue)) 
        [queue,tmp0] = dequeue(queue);
        neighbor = enqueue(neighbor, tmp0);
        
        for s = 1:size(move,1)
            tmp1 = tmp0 + move(s,:);
            if(tmp1(1)>=iMin && tmp1(1)<=iMax && ...
               tmp1(2)>=jMin && tmp1(2)<=jMax && ...
               tmp1(3)>=kMin && tmp1(3)<=kMax && ...
               ~isVisit(tmp1(1),tmp1(2),tmp1(3)))
           
                queueNext = enqueue(queueNext, tmp1);
                isVisit(tmp1(1),tmp1(2),tmp1(3)) = true;
            end
        end
    end
    queue = queueNext;
end

neighbor = [neighbor; queue];

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Queue operation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function queue = enqueue(queue, x)

isDup = ismember(x,queue,'rows');
if(~isDup)
    N = size(queue,1);
    queue(N+1,:) = x;
end

end


function [queue,x] = dequeue(queue)

if(size(queue,1)==0)
    error('Empty queue.');
end
x = queue(1,:);
queue(1,:) = [];

end

% END OF FILE