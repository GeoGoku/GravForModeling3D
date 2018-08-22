function [meshInfo, cell_n, cell_e, cell_v, cell_dn, cell_de, cell_dv] = Mesh3D(meshFile)
%MESH3D transfers the meshFile into a 3d case.
%   MESH3D transfers the geometric parameters of 3D geological body from UBC-MeshTools3d format to xyz format
%   based on meshFile.
%
%   More info please refer to UBC-MeshTools3d manual.
%
%   Input:
%             'meshFile':                    filename of mesh
%
%   Output:
%             'cell_n':                         center coordinates of cells in x/north direction
%             'cell_e':                         center coordinates of cells in y/east direction
%             'cell_v':                         center coordinates of cells in z/down direction
%             'cell_dn':                       size of cells in x/north direction
%             'cell_de':                       size of cells in y/east direction
%             'cell_dv':                       size of cells in z/down direction
%             
%   Example:
%     meshFile = 'cubic';
%     [meshInfo, cell_n, cell_e, cell_v, cell_dn, cell_de, cell_dv] = Mesh3D(meshFile)
%
%   T. Chen 20-May-2016.
%   Upgrade in 15-Mar-2017.
%   Copyright 2016-2017 T. Chen.

% Check arguments.
nbIn = nargin;
switch nbIn
    case {1}
    case {0}
        error('Mesh3D: FunctionInput: NotEnough_ArgNum');
    otherwise
        error('Mesh3D: FunctionInput: TooMany_ArgNum');
end

[fileID,errmsg] = fopen(meshFile,'r');

if fileID < 0
    disp(['Mesh3D: FileInput_', errmsg]);
end

a1 = fgets(fileID);                                                    % 剖分个数(NE,NN,NZ)
a2 = fgets(fileID);                                                    % 西顶角坐标(E0,N0,Z0)
a3 = fgets(fileID);                                                    % 东向剖分信息
a4 = fgets(fileID);                                                    % 北向剖分信息
a5 = fgets(fileID);                                                    % 垂向剖分信息
fclose(fileID);

% num of mesh cells in east, north and down direction
% 注意：mesh文件中的方向为东向、北向、垂直方向，具体可参考UBC系列软件说明书
b = regexp(a1,'\d*\.?\d*','match');
NE = str2num(cell2mat(b(1)));
NN = str2num(cell2mat(b(2)));
NV = str2num(cell2mat(b(3)));
meshInfo = [NN NE NV];

% 西南顶角坐标（top-south-west）
b = regexp(a2,'(-\d|\d*)*','match');
min_e = str2num(cell2mat(b(1)));
min_n = str2num(cell2mat(b(2)));
min_v = str2num(cell2mat(b(3)));

% 东向-方向剖分信息（east）
b = regexp(a3,'\d*\.?\d*','match');
mesh_ne = zeros(1,size(b,2)/2);       % num of cells in east direction/东向-方向网格剖分个数
mesh_de = zeros(1,size(b,2)/2);       % size of cells in east direction/东向-方向网格剖分间距
for i=1:size(b,2)/2
    mesh_ne(i) = str2num(cell2mat(b(2*(i-1)+1)));
    mesh_de(i) = str2num(cell2mat(b(2*(i-1)+2)));
end

% 北向-方向剖分信息（north）
b = regexp(a4,'\d*\.?\d*','match');
mesh_nn = zeros(1,size(b,2)/2);       % num of cells in north direction/北向-方向网格剖分个数
mesh_dn = zeros(1,size(b,2)/2);       % size of cells in north direction/北向-方向网格剖分间距
for i=1:size(b,2)/2
    mesh_nn(i) = str2num(cell2mat(b(2*(i-1)+1)));
    mesh_dn(i) = str2num(cell2mat(b(2*(i-1)+2)));
end

% 垂向-方向剖分（vertical/depth）
b = regexp(a5,'\d*\.?\d*','match');
mesh_nv = zeros(1,size(b,2)/2);       % num of cells in down direction/垂向-方向网格剖分个数
mesh_dv = zeros(1,size(b,2)/2);       % size of cells in down direction/垂向-方向网格剖分间距
for i=1:size(b,2)/2
    mesh_nv(i) = str2num(cell2mat(b(2*(i-1)+1)));
    mesh_dv(i) = str2num(cell2mat(b(2*(i-1)+2)));
end

% east
cell_de = [];
for i = 1 : size(mesh_ne, 2)
    add = mesh_de(i) * ones(1, mesh_ne(i));
    cell_de = [cell_de add];
end

cell_e = cell_de;
for i = 1 : length(cell_e)
    cell_e(i) = min_e + sum(cell_de(1:i)) - cell_de(i) * 0.5;
end

% north
cell_dn = [];
for i = 1 : size(mesh_nn, 2)
    add = mesh_dn(i) * ones(1, mesh_nn(i));
    cell_dn = [cell_dn add];
end

cell_n = cell_dn;
for i = 1 : length(cell_n)
    cell_n(i) = min_n + sum(cell_dn(1:i)) - cell_dn(i) * 0.5;
end

% vertical
cell_dv = [];
for i = 1 : size(mesh_nv, 2)
    add = mesh_dv(i) * ones(1, mesh_nv(i));
    cell_dv = [cell_dv add];
end

cell_v = cell_dv;
for i = 1 : length(cell_v)
    cell_v(i) = min_v + sum(cell_dv(1:i)) - cell_dv(i) * 0.5;
end

% The model ordering is performed first in the z/vertical direction (top-to-bottom),
% then in the easting, and finally in the northing.
cell_e = repmat(cell_e, NV, 1, NN);
cell_n = reshape(cell_n, 1, 1, NN);
cell_n = repmat(cell_n, NV, NE, 1);
cell_v = repmat(cell_v', 1, NE, NN);
cell_de = repmat(cell_de, NV, 1, NN);
cell_dn = reshape(cell_dn, 1, 1, NN);
cell_dn = repmat(cell_dn, NV, NE, 1);
cell_dv = repmat(cell_dv', 1, NE, NN);

cell_e = reshape(cell_e, NV * NE * NN, 1);
cell_n = reshape(cell_n, NV * NE * NN, 1);
cell_v = reshape(cell_v, NV * NE * NN, 1);
cell_de = reshape(cell_de, NV * NE * NN, 1);
cell_dn = reshape(cell_dn, NV * NE * NN, 1);
cell_dv = reshape(cell_dv, NV * NE * NN, 1);

end