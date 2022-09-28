function [strand,base2node] = transformMat_original(dnaTop,strand,HJE,CanDoPath,filename)


%% Load data
%deformTrans = cell(nHelix,1);  % deformed translation
load([CanDoPath, '\', filename, '_FEModel.mat']);
load([CanDoPath, '\', filename, '_Displ.mat']);
deformedCoordinates = FEModel.node(:,1:3) + displ(:,1:3);
clear FEModel displ
initialRotations = calcInitialRotations(CanDoPath,filename);
deformedRotations = calcDeformedRotations(CanDoPath,filename);


% Change the sequence of the node from ADINA convention to Tiamat topology
% base ID convention
nStrand = numel(strand);
nHJE = numel(HJE);
nBase = numel(dnaTop);
base2node = zeros(nBase,1);

iNode = 0;
for i = 1:nHJE
    for j = 1:numel(HJE(i).idSeq_main)
        iNode = iNode + 1;
        base2node(HJE(i).idSeq_main(j)) = iNode;
    end
end
iNode = 0;
for i = 1:nHJE
    for j = 1:numel(HJE(i).idSeq_comp)
        iNode = iNode + 1;
        base2node(HJE(i).idSeq_comp(j)) = iNode;
    end
end
%assert(isempty(find(base2node<=0, 1)));
for i = 1:numel(base2node)
    if(base2node(i)==0 && dnaTop(i).across>=0)
        error('Undefined base-pair position and orientation.');
    end
end

for i = 1:nStrand
    for j = 1:numel(strand(i).tour)
        strand(i).R{j} = [];
        strand(i).d{j} = [];
        iNode = base2node(strand(i).tour(j));
        if(iNode > 0)
            R0 = initialRotations(:,:,iNode);
            R1 = deformedRotations(:,:,iNode);
            d = deformedCoordinates(iNode,:);
            strand(i).R{j} = R1*R0;
            strand(i).d{j} = d';
        end
    end
end

% s2 = 0;
% for i = 1:nHelix
%     seqInd = deformList.helixSeq(i);
%     if(seqInd >= 1)
%         s1 = s2+1;
%         s2 = s2+lenHelix(seqInd);
%         initialRot{seqInd} = initialRotations(:,:,s1:s2);
%         deformRot{seqInd} = deformedRotations(:,:,s1:s2);
%         deformTrans{seqInd} = deformedCoordinates(s1:s2,:);
%     end
% end
% clear initialCoordinates deformedCoordinates deformedRotations;


%% Calculate the rotation/translation
% [scafTourLineR,scafTourLineD] = calcGeometry(deformList.scafTourLine,initialRot,deformRot,deformTrans);
% [scafTourLoopR,scafTourLoopD] = calcGeometry(deformList.scafTourLoop,initialRot,deformRot,deformTrans);
% [stapTourLineR,stapTourLineD] = calcGeometry(deformList.stapTourLine,initialRot,deformRot,deformTrans);
% [stapTourLoopR,stapTourLoopD] = calcGeometry(deformList.stapTourLoop,initialRot,deformRot,deformTrans);

end


% Calculate rotation matrix R and translation vector D
function [finalR,finalD] = calcGeometry(tour,initialRot,deformRot,deformTrans)

nTour = length(tour);
finalR = cell(nTour,1);
finalD = cell(nTour,1);

for i = 1:nTour
    lenTour = size(tour{i},1);
    finalR{i} = zeros(3,3,lenTour);
    finalD{i} = zeros(3,1,lenTour);
    for j = 1:lenTour
        % Which helix?
        h = tour{i}(j,1);
        % Which position?
        p = tour{i}(j,2);
        % Initial rotation
        R0 = initialRot{h}(:,:,p);
        % Deformed rotation
        R1 = deformRot{h}(:,:,p);
        
        % Deformed translation
        D = deformTrans{h}(p,:);
        
        % Save results
        finalR{i}(:,:,j) = R1*R0;
        finalD{i}(:,:,j) = D';
        
        % Which tour for staple?
        %i1 = deformList.stapTourMark{h}(p,1);
        % Which position in the tour?
        %j1 = deformList.stapTourMark{h}(p,2);
        
        % Save results
        %if(i1>0)
        %    finalR{i1}(:,:,j1) = R1*R0;
        %    finalD{i1}(:,:,j1) = D';
        %else
        %    finalR{-i1}(:,:,j1) = R1*R0;
        %    finalD{-i1}(:,:,j1) = D';
        %end
        
        %fprintf('Circular scaffold %d-%d in structure %d-%d, initial rotation=%.2f, deformed rotation= %.2f %.2f %.2f\n', ...
        %        i, j, h, p, initialRot{h}(p), deformRot{h}(p,1)/pi*180, deformRot{h}(p,2)/pi*180, deformRot{h}(p,3)/pi*180);
    end
end

end


% Rotation about X-axis
function R = Rx(deg)
rad = deg/180*pi;
R = [1 0 0; 0 cos(rad) -sin(rad); 0 sin(rad) cos(rad)];
end

% Rotation about Y-axis
function R = Ry(deg)
rad = deg/180*pi;
R = [cos(rad) 0 sin(rad); 0 1 0; -sin(rad) 0 cos(rad)];
end

% Rotation about Z-axis
function R = Rz(deg)
rad = deg/180*pi;
R = [cos(rad) -sin(rad) 0; sin(rad) cos(rad) 0; 0 0 1];
end


function M = rotMat(theta)
% ----- Input -----
% theta is a 1-by-3 row vector, which is one row in CanDo output
% deformedRotation.txt
% The magnitude of theta (phi in the program) is the rotation angle
% The direction of theta (u in the program) is the rotation axis
% ----- Output -----
% M is a 3-by-3 rotation matrix.
% ----- Example -----
% Let x be a 3-by-1 column vector, now rotate it by theta:
% x = M * x
%

phi = norm(theta);  % phi is the rotation angle around the rotation axis
if(0==phi)
    M = eye(3);
    return;
end

u = theta/phi;   % u is the rotation axis

% Create the rotation matrix, see
% http://en.wikipedia.org/wiki/Rotation_matrix for the formula
M = zeros(3);

M(1,1) = cos(phi) + (1-cos(phi))*u(1)^2;
M(1,2) = (1-cos(phi))*u(1)*u(2) - u(3)*sin(phi);
M(1,3) = (1-cos(phi))*u(1)*u(3) + u(2)*sin(phi);

M(2,1) = (1-cos(phi))*u(1)*u(2) + u(3)*sin(phi);
M(2,2) = cos(phi) + (1-cos(phi))*u(2)^2;
M(2,3) = (1-cos(phi))*u(2)*u(3) - u(1)*sin(phi);

M(3,1) = (1-cos(phi))*u(1)*u(3) - u(2)*sin(phi);
M(3,2) = (1-cos(phi))*u(2)*u(3) + u(1)*sin(phi);
M(3,3) = cos(phi) + (1-cos(phi))*u(3)^2;

end


function M = rotMat2(theta)
% Test case: 
% jsonPath = 'C:\ADINA_run\MURI';
% filename = '6-helix-bend';

% Convert radian to angle
a1 = theta(1)*180/pi;
a2 = theta(2)*180/pi;
a3 = theta(3)*180/pi;

% Euler angles
%M = Rz(a1)*Ry(a2)*Rz(a3);  % Very wrong for a 6-Helix bundle bending about
% z-axis with a negative angle
% RPY
%M = Rz(a3)*Rx(a1)*Ry(a2);

M = Rz(a3)*Ry(a2)*Rx(a1);
%M = Rx(a1)*Ry(a2)*Rz(a3);   % Very wrong for a 6-Helix bundle bending about
% z-axis with a negative angle
end


function R = calcInitialRotations(CanDoPath,filename)

load(fullfile(CanDoPath, [filename, '_FEModel.mat']));
nbp = size(FEModel.node,1);
R = zeros(3,3,nbp);

for i = 1:nbp
    e1 = FEModel.nodalTriad(i,:,1);
    e1 = e1/norm(e1);
    e2 = FEModel.nodalTriad(i,:,2);
    e2 = e2/norm(e2);
    e3 = cross(e1,e2);
    R(:,:,i) = [e3', e1', e2'];
end

end


function R = calcDeformedRotations(CanDoPath,filename)
maxNQuat = 1e7;
% Read quaternions
triadID = zeros(maxNQuat,1);
quat0 = [];
quat1 = zeros(maxNQuat,4);
quat2 = zeros(maxNQuat,4);

fid = fopen([CanDoPath, '\', filename, '.por']);
tline = fgetl(fid);

ind = 0;
while ischar(tline)
    if(strcmp(strtrim(tline),'CURTRIAD'))
        ind = ind+1;
        if(ind>maxNQuat)
            error('The number of quaternions overflow.');
        end
        
        tline = fgetl(fid);
        triadID(ind) = str2double(tline);
        tline = fgetl(fid);
        quat1(ind,:) = [str2double(tline(1:20)) str2double(tline(21:40)) str2double(tline(41:60)) str2double(tline(61:80))];
        tline = fgetl(fid);
        quat2(ind,:) = [str2double(tline(1:20)) str2double(tline(21:40)) str2double(tline(41:60)) str2double(tline(61:80))];
    elseif(strcmp(strtrim(tline),'INITRIAD'))
        tline = fgetl(fid);
        quat0 = [quat0; str2double(tline(1:20)) str2double(tline(21:40)) str2double(tline(41:60)) str2double(tline(61:80))];
    else
        tline = fgetl(fid);
    end
end
fclose(fid);
triadID(ind+1:maxNQuat) = [];
quat1(ind+1:maxNQuat,:) = [];
quat2(ind+1:maxNQuat,:) = [];

if(~isempty(find(isnan(triadID),1)) || ~isempty(find(isinf(triadID),1)))
    error('Illegal value in triadID');
end
if(~isempty(find(isnan(quat1),1)) || ~isempty(find(isinf(quat1),1)))
    error('Illegal value in quat1');
end
if(~isempty(find(isnan(quat2),1)) || ~isempty(find(isinf(quat2),1)))
    error('Illegal value in quat2');
end  


load([CanDoPath, '\', filename, '_FEModel.mat']);
% Discussion with Do-Nyun on 08/14/2012:
% FEModel.beamsNormal and FEModel.beamsNick store element information.
% Columns 1-3: ElementID, FirstNodeID, SecondNodeID
% The number of rows of FEModel.node is the number of base pairs.
% See the email forwarded by Do-Nyun on 08/14/2012
beamsElement = FEModel.beamsNormal;
if(isfield(FEModel, 'beamsNick'))
    beamsElement = [beamsElement; FEModel.beamsNick];
end
nbp = size(FEModel.node,1);
nElement = size(beamsElement,1);
if(isfield(FEModel, 'rigidLinks'))
    nRigidLink = size(FEModel.rigidLinks,1);
else
    nRigidLink = 0;
end

if(length(unique([beamsElement(:,2);beamsElement(:,3)])) ~= length(FEModel.node))
    error('Some nodes are exclusively connected with elements that are neither normal beam or nick.');
end
[~,indLast] = ismember(beamsElement(:,1),triadID, 'legacy');
if(min(indLast) <= 0)
    error('Empty beam ID.');
end
triadID = triadID(indLast);
if(length(triadID) ~= length(unique(triadID)))
    error('Duplication in triadID');
end
quat1 = quat1(indLast,:);
quat2 = quat2(indLast,:);

%if(length(triadID) ~= length(unique(triadID)))
%    % Duplication in triadID
%    tmp = length(triadID)-nElement-nRigidLink+1 : length(triadID);
%    triadID = triadID(tmp);
%    quat1 = quat1(tmp,:);
%    quat2 = quat2(tmp,:);
%    if(length(triadID) ~= length(unique(triadID)))
%        error('Duplication in triadID');
%    end
%end
if(~issorted(triadID))  
    %error('Indices of the triads are not sorted.');
    % If FEModel.beamsNick is not empty, triadID is not sorted because the triad is [FEModel.beamsNormal;FEModel.beamsNick;FEModel.rigidLinks]
    [triadID,si] = sort(triadID);
    quat0 = quat0(si,:);
    quat1 = quat1(si,:);
    quat2 = quat2(si,:);
end


quatNode = zeros(nbp,4);
quatFilled = zeros(nbp,1);
quatNodeInit = zeros(nbp,4);
quatFilledInit = zeros(nbp,1);
for i = 1:size(beamsElement,1)
    elementID = beamsElement(i,1);
    node1ID = beamsElement(i,2);
    node2ID = beamsElement(i,3);
    
    if(1==quatFilled(node1ID) && norm(quatNode(node1ID,:)-quat1(elementID,:))>1e-4)
        error('Conflicting quaternions after deformation');
    end
    quatNode(node1ID,:) = quat1(elementID,:);
    quatFilled(node1ID) = 1;
    
    if(1==quatFilled(node2ID) && norm(quatNode(node2ID,:)-quat2(elementID,:))>1e-4)
        error('Conflicting quaternions after deformation');
    end
    quatNode(node2ID,:) = quat2(elementID,:);
    quatFilled(node2ID) = 1;
    
    if(1==quatFilledInit(node1ID) && norm(quatNodeInit(node1ID,:)-quat0(elementID,:))>1e-4)
        error('Conflicting quaternions before deformation');
    end
    quatNodeInit(node1ID,:) = quat0(elementID,:);
    quatFilledInit(node1ID) = 1;
    
    if(1==quatFilledInit(node2ID) && norm(quatNodeInit(node2ID,:)-quat0(elementID,:))>1e-4)
        error('Conflicting quaternions before deformation');
    end
    quatNodeInit(node2ID,:) = quat0(elementID,:);
    quatFilledInit(node2ID) = 1;
end
if(min(quatFilled)<1 || min(quatFilledInit)<1)
    error('Some nodes are not assigned.');
end

% Convert quaternions to DCMs
DCMinit = quat2dcm(quatNodeInit);
DCMcurr = quat2dcm(quatNode);
R = zeros(3,3,nbp);
for i = 1:nbp
    R(:,:,i) = DCMcurr(:,:,i) * DCMinit(:,:,i)';
end

end