function [Mesh_Cellgroup] = CellMerge1D(mesh1, mesh2, mesh3, mesh4, mesh5, mesh6, model, threshold)
% CELLMERGE1D performs 1d cell mergence
%   CELLMERGE1D performs 1d cell mergence using
%   mesh1, mesh2, mesh3, mesh4, mesh5, mesh6, model, threshold
%
%  Input：    
%         mesh1                 待融合的网格单元—x/north coordinate
%         mesh2                 待融合的网格单元—y/east coordinate
%         mesh3                 待融合的网格单元—z/down coordinate
%         mesh4                 待融合的网格单元—size in x/north direction
%         mesh5                 待融合的网格单元—size in y/east direction
%         mesh6                 待融合的网格单元—size in z/down direction
%         model                  待融合的网格单元—model
%         threshold             网格单元融合阈值
%
%  Output：
%         mesh_cellgroup    merged cell groups/融合后的网格单元组
%
%   T. Chen 24-Sep-2016.
%   Upgrade in 23-Nov-2016.
%   Copyright 2016 T. Chen.
%   Notes: The CELLMERGE1D requires the mesh composed of rectargular cells (they could be unequal interval).
%   Notes: CELLMERGE1D 以列融合为例进行编程计算。以行融合的变换参数获取。 

% Convert vector to column vector/将向量转化为列向量
mesh1  = mesh1(:)';  mesh2  = mesh2(:)';
mesh3  = mesh3(:)';  mesh4  = mesh4(:)';
mesh5  = mesh5(:)';  mesh6  = mesh6(:)';
model  = model(:)';

% Remove non-zero values/根据密度值将非零值去掉
Index_Nonzeros = find(model ~= 0);

% Find the backward difference of the density vector/对密度向量求后向差分
B = [diff(model) model(end)];
% B(abs(B) <= threshold) = 0;

% Find the forward difference of the density vector/对密度向量求前向差分
C = fliplr(model);
D = [diff(C) C(end)];
% D(abs(D) <= threshold) = 0;
D = fliplr(D);

% calculate the flag value of cells and its location/计算网格单元标志值（标志值为0）及其位置
Index_Position = find(B.*D == 0);
Mesh3_Cellgroup = model(Index_Position);
Index_Cellgroup = Index_Position;
if ismember(0, Mesh3_Cellgroup)
    Index_Position = find(Mesh3_Cellgroup ~= 0);
    Mesh3_Cellgroup = Mesh3_Cellgroup(Index_Position);
    Index_Cellgroup = Index_Cellgroup(Index_Position);
end

% size of cellgroup
Res_Cell = setdiff(Index_Nonzeros, Index_Cellgroup);
j = 1;

if ~isempty(Mesh3_Cellgroup)
    % find location of cell group/找出网格单元组位置
    E = gradient(Mesh3_Cellgroup);
    F = find(E);
    F = [1 F length(E)];

    B = [diff(Index_Cellgroup) Index_Cellgroup(end)];
    FB = find(B == 1);
    C = fliplr(Index_Cellgroup);
    D = [diff(C) C(end)];
    D = fliplr(D);
    FD = find(D == -1);
    BuD = intersect(FB, FD);
    BnD = [FB FD];
    G = setdiff(BnD,BuD);

    index = unique([F G]);

    isgpuArray = whos('mesh1');
    if strcmp(isgpuArray.class,'gpuArray')
        Mesh_Cellgroup = gpuArray.zeros(length(Res_Cell) + length(index)/2, 7);
    else
        Mesh_Cellgroup = zeros(length(Res_Cell) + length(index)/2, 7);
    end

    % cell merge
    for i = 1:length(index)/2
        Index_Merge = Index_Cellgroup(index(2 * i - 1):index(2 * i));
        Mesh_Cellgroup(j, 1) = mesh1(Index_Merge(1));                                                                                          % geometric parameter merge/x/north
        Mesh_Cellgroup(j, 2) = mesh2(Index_Merge(1));                                                                                          % geometric parameter merge/y/east
        Mesh_Cellgroup(j, 3) = mesh3(Index_Merge(1)) + sum(mesh6(Index_Merge(2):Index_Merge(end)))/2;         % geometric parameter merge/z/down

        Mesh_Cellgroup(j, 4) = mesh4(Index_Merge(1));                                                                                          % geometric parameter merge/dx/north
        Mesh_Cellgroup(j, 5) = mesh5(Index_Merge(1));                                                                                          % geometric parameter merge/dy/east
        Mesh_Cellgroup(j, 6) = sum(mesh6(Index_Merge));                                                                                     % geometric parameter merge/dz/down
        
        Mesh_Cellgroup(j, 7) = sum(model(Index_Merge))/length(Index_Merge);                                                    % physical parameter merge/density
        j = j + 1;
    end
else
    isgpuArray = whos('mesh1');
    if strcmp(isgpuArray.class,'gpuArray')
        Mesh_Cellgroup = gpuArray.zeros(length(Res_Cell), 7);
    else
        Mesh_Cellgroup = zeros(length(Res_Cell), 7);
    end
end

% estaiblish cellgroups for each unmerged cell/将每一个没有融合的网格单元单独建立一个单元组
for i = 1 : length(Res_Cell)
    Mesh_Cellgroup(j, 1) = mesh1(Res_Cell(i));
    Mesh_Cellgroup(j, 2) = mesh2(Res_Cell(i));
    Mesh_Cellgroup(j, 3) = mesh3(Res_Cell(i));
    Mesh_Cellgroup(j, 4) = mesh4(Res_Cell(i));
    Mesh_Cellgroup(j, 5) = mesh5(Res_Cell(i));
    Mesh_Cellgroup(j, 6) = mesh6(Res_Cell(i));
    Mesh_Cellgroup(j, 7) = model(Res_Cell(i));
    j = j + 1;
end

end