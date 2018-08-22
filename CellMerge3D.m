function [cellgroup_n, cellgroup_e, cellgroup_v, cellgroup_dn, cellgroup_de, cellgroup_dv, cellgroup_model] = CellMerge3D(meshFile, modelFile, threshold, parFlag)
%CELLMERGE3D performs 3d cell mergence
%   CELLMERGE3D performs 3d cell mergence using specified parameters
%
%   Input：    
%          meshFile                                        file name of mesh
%          modelFile                                      file name of model
%          threshold                                      threshold for cell merging
%          parFlag                                          flag of parallel computing (1 serial, 2 for CPU parallel, 3 for GPU parallel）
%
%   Output：
%          cellgroup_n                                   center coordinates of cellgroups in x/north direction
%          cellgroup_e                                   center coordinates of cellgroups in y/east direction
%          cellgroup_v                                   center coordinates of cellgroups in z/down direction
%          cellgroup_dn                                 size of cellgroups in x/north direction
%          cellgroup_de                                 size of cellgroups in y/east direction
%          cellgroup_dv                                 size of cellgroups in z/down direction
%          cellgroup_model                           physical parameters of cellgroups
%
%   T. Chen 27-Sep-2016.
%   Upgrade in 16-Feb-2017.
%   Copyright 2016-2017 T. Chen.

% read model and mesh info
[meshInfo, mesh_x_north, mesh_y_east, mesh_z_down, mesh_dx_north, mesh_dy_east, mesh_dz_down] = Mesh3D(meshFile);
mesh_model = importdata(modelFile);
modelsize1 = length(nonzeros(mesh_model));

% reshape the mesh string/对原始的模型网格单元数据进行三维矩阵变形
mesh_x_north  = reshape(mesh_x_north, meshInfo(3), meshInfo(2), meshInfo(1));  % geometric matrix (x/north coordinate)
mesh_y_east   = reshape(mesh_y_east, meshInfo(3), meshInfo(2), meshInfo(1));   % geometric matrix (y/east coordinate)
mesh_z_down   = reshape(mesh_z_down, meshInfo(3), meshInfo(2), meshInfo(1));   % geometric matrix (z/down coordinate)
mesh_model    = reshape(mesh_model, meshInfo(3), meshInfo(2), meshInfo(1));    % physical parameter matrix (model)
mesh_dx_north = reshape(mesh_dx_north, meshInfo(3), meshInfo(2), meshInfo(1)); % geometric matrix (size in x/north direction)
mesh_dy_east  = reshape(mesh_dy_east, meshInfo(3), meshInfo(2), meshInfo(1));  % geometric matrix (size in y/east direction)
mesh_dz_down  = reshape(mesh_dz_down, meshInfo(3), meshInfo(2), meshInfo(1));  % geometric matrix (size in z/down direction)

% 统计模型在北、东、垂三个方向的尺度大小。注意矩阵索引和实际地质体单元的差异
[rows_down, cols_east, pags_north] = ind2sub(size(mesh_model), find(mesh_model)); % 找出非零值，减少计算量
sacle_north = unique(pags_north);
sacle_east = unique(cols_east);
sacle_down = unique(rows_down);
% 将矩阵分解为一维的向量然后网格合并
% 按照尺度来判单元合并方向
% (1) 北向尺度大，北向合并
% (2) 东向尺度大，东向合并
% (3) 垂向尺度大，垂向合并
% (4) 尺度相同，按照垂向合并，因为MATLAB中高维矩阵存储顺序为：列优先，有利于增速。
% merge in east direction/东向合并    
if length(sacle_east)>length(sacle_down) && length(sacle_east)>length(sacle_north)
    index = [pags_north, rows_down];
    index = unique(index, 'rows');
    ct = size(index,1);
    mesh_cellgroup = cell(ct,1);
    switch parFlag
        case 1
            for i = 1 : ct            % attention: row-down; col-east; page-north, this is determined by the storage order
                mesh1  = mesh_x_north(index(i, 2), :, index(i, 1));
                mesh2  = mesh_y_east(index(i, 2), :, index(i, 1));
                mesh3  = mesh_z_down(index(i, 2), :, index(i, 1));
                mesh4  = mesh_dx_north(index(i, 2), :, index(i, 1));
                mesh5  = mesh_dy_east(index(i, 2), :, index(i, 1));
                mesh6  = mesh_dz_down(index(i, 2), :, index(i, 1));
                model1  =  mesh_model(index(i, 2), :, index(i, 1));
                % The order of the parameters changes to fit cellmerge1d. Because cellmerge1d is written in a vertical merge.
                mesh_cellgroup{i} = CellMerge1D(mesh1(:), mesh3(:), mesh2(:), mesh4(:), mesh6(:), mesh5(:), model1(:), threshold);
            end
        case 2
            parfor i = 1 : ct
                mesh1  = mesh_x_north(index(i, 2), :, index(i, 1));
                mesh2  = mesh_y_east(index(i, 2), :, index(i, 1));
                mesh3  = mesh_z_down(index(i, 2), :, index(i, 1));
                mesh4  = mesh_dx_north(index(i, 2), :, index(i, 1));
                mesh5  = mesh_dy_east(index(i, 2), :, index(i, 1));
                mesh6  = mesh_dz_down(index(i, 2), :, index(i, 1));
                model1  =  mesh_model(index(i, 2), :, index(i, 1));
                % The order of the parameters changes to fit cellmerge1d. Because cellmerge1d is written in a vertical merge.
                mesh_cellgroup{i} = CellMerge1D(mesh1(:), mesh3(:), mesh2(:), mesh4(:), mesh6(:), mesh5(:), model1(:), threshold);
            end
        case 3
            mesh_y_east = gpuArray(mesh_y_east); mesh_z_down = gpuArray(mesh_z_down);
            mesh_model = gpuArray(mesh_model); mesh_dy_east = gpuArray(mesh_dy_east);
            mesh_dz_down = gpuArray(mesh_dz_down);
            for i = 1 : ct
                mesh1  = mesh_x_north(index(i, 2), :, index(i, 1));
                mesh2  = mesh_y_east(index(i, 2), :, index(i, 1));
                mesh3  = mesh_z_down(index(i, 2), :, index(i, 1));
                mesh4  = mesh_dx_north(index(i, 2), :, index(i, 1));
                mesh5  = mesh_dy_east(index(i, 2), :, index(i, 1));
                mesh6  = mesh_dz_down(index(i, 2), :, index(i, 1));
                model1  =  mesh_model(index(i, 2), :, index(i, 1));
                % The order of the parameters changes to fit cellmerge1d. Because cellmerge1d is written in a vertical merge.
                mesh_cellgroup{i} = CellMerge1D(mesh1(:), mesh3(:), mesh2(:), mesh4(:), mesh6(:), mesh5(:), model1(:), threshold);
            end
    end
    mesh_cellgroup = cell2mat(mesh_cellgroup);
    cellgroup_model = mesh_cellgroup(:, 7);
    cellgroup_n = mesh_cellgroup(:, 1);
    cellgroup_e = mesh_cellgroup(:, 3);
    cellgroup_v = mesh_cellgroup(:, 2);
    cellgroup_dn = mesh_cellgroup(:, 4);
    cellgroup_de = mesh_cellgroup(:, 6);
    cellgroup_dv = mesh_cellgroup(:, 5);
% merge in north direction/北向合并
elseif length(sacle_north)>length(sacle_down) && length(sacle_north)>length(sacle_east)
    index = [cols_east, rows_down];
    index = unique(index, 'rows');
    ct = size(index,1);
    mesh_cellgroup = cell(ct,1);
    switch parFlag
        case 1
            for i = 1 : ct            % attention: row-down; col-east; page-north, this is determined by the storage order
                % The order of the parameters changes to fit cellmerge1d. Because cellmerge1d is written in a vertical merge.
                mesh_cellgroup{i} = CellMerge1D(mesh_z_down(index(i, 2),  index(i, 1), :), mesh_y_east(index(i, 2),  index(i, 1), :), ...
                    mesh_x_north(index(i, 2),  index(i, 1), :), mesh_dz_down(index(i, 2),  index(i, 1), :), ...
                    mesh_dy_east(index(i, 2),  index(i, 1), :), mesh_dx_north(index(i, 2),  index(i, 1), :), ...
                    mesh_model(index(i, 2),  index(i, 1), :), threshold);
            end
        case 2
            parfor i = 1 : ct
                mesh_cellgroup{i} = CellMerge1D(mesh_z_down(index(i, 2),  index(i, 1), :), mesh_y_east(index(i, 2),  index(i, 1), :), ...
                    mesh_x_north(index(i, 2),  index(i, 1), :), mesh_dz_down(index(i, 2),  index(i, 1), :), ...
                    mesh_dy_east(index(i, 2),  index(i, 1), :), mesh_dx_north(index(i, 2),  index(i, 1), :), ...
                    mesh_model(index(i, 2),  index(i, 1), :), threshold);
            end
        case 3
            mesh_x_north = gpuArray(mesh_x_north);mesh_y_east = gpuArray(mesh_y_east); mesh_z_down = gpuArray(mesh_z_down);
            mesh_model = gpuArray(mesh_model); mesh_dx_north = gpuArray(mesh_dx_north); mesh_dy_east = gpuArray(mesh_dy_east);
            mesh_dz_down = gpuArray(mesh_dz_down);
            for i = 1 : ct
                mesh_cellgroup{i} = CellMerge1D(mesh_z_down(index(i, 2),  index(i, 1), :), mesh_y_east(index(i, 2),  index(i, 1), :), ...
                    mesh_x_north(index(i, 2),  index(i, 1), :), mesh_dz_down(index(i, 2),  index(i, 1), :), ...
                    mesh_dy_east(index(i, 2),  index(i, 1), :), mesh_dx_north(index(i, 2),  index(i, 1), :), ...
                    mesh_model(index(i, 2),  index(i, 1), :), threshold);
            end
    end
    mesh_cellgroup = cell2mat(mesh_cellgroup);
    cellgroup_model = mesh_cellgroup(:, 7);
    cellgroup_n = mesh_cellgroup(:, 3);
    cellgroup_e = mesh_cellgroup(:, 2);
    cellgroup_v = mesh_cellgroup(:, 1);
    cellgroup_dn = mesh_cellgroup(:, 6);
    cellgroup_de = mesh_cellgroup(:, 5);
    cellgroup_dv = mesh_cellgroup(:, 4);
else  %length(sacle_down)>length(sacle_east) && length(sacle_down)>length(sacle_north)
    index = [pags_north, cols_east];
    index = unique(index, 'rows');
    ct = size(index,1);
    mesh_cellgroup = cell(ct,1);
    switch parFlag
        case 1
            for i = 1 : ct           % attention: row-down; col-east; page-north, this is determined by the storage order
                mesh_cellgroup{i} = CellMerge1D(mesh_x_north(:, index(i, 2), index(i, 1)), ...
                    mesh_y_east(:, index(i, 2), index(i, 1)), mesh_z_down(:, index(i, 2), index(i, 1)),...
                    mesh_dx_north(:, index(i, 2), index(i, 1)), mesh_dy_east(:, index(i, 2), index(i, 1)), ...
                    mesh_dz_down(:, index(i, 2), index(i, 1)), mesh_model(:, index(i, 2), index(i, 1)), threshold);
            end
        case 2
            parfor i = 1 : ct
                mesh_cellgroup{i} = CellMerge1D(mesh_x_north(:, index(i, 2), index(i, 1)), ...
                    mesh_y_east(:, index(i, 2), index(i, 1)), mesh_z_down(:, index(i, 2), index(i, 1)),...
                    mesh_dx_north(:, index(i, 2), index(i, 1)), mesh_dy_east(:, index(i, 2), index(i, 1)), ...
                    mesh_dz_down(:, index(i, 2), index(i, 1)), mesh_model(:, index(i, 2), index(i, 1)), threshold);
            end
        case 3
            mesh_x_north = gpuArray(mesh_x_north);mesh_y_east = gpuArray(mesh_y_east); mesh_z_down = gpuArray(mesh_z_down);
            mesh_model = gpuArray(mesh_model); mesh_dx_north = gpuArray(mesh_dx_north); mesh_dy_east = gpuArray(mesh_dy_east);
            mesh_dz_down = gpuArray(mesh_dz_down);
            for i = 1 : ct
                mesh_cellgroup{i} = CellMerge1D(mesh_x_north(:, index(i, 2), index(i, 1)), ...
                    mesh_y_east(:, index(i, 2), index(i, 1)), mesh_z_down(:, index(i, 2), index(i, 1)),...
                    mesh_dx_north(:, index(i, 2), index(i, 1)), mesh_dy_east(:, index(i, 2), index(i, 1)), ...
                    mesh_dz_down(:, index(i, 2), index(i, 1)), mesh_model(:, index(i, 2), index(i, 1)), threshold);
            end
    end
    mesh_cellgroup = cell2mat(mesh_cellgroup);
    cellgroup_model = mesh_cellgroup(:, 7);
    cellgroup_n = mesh_cellgroup(:, 1);
    cellgroup_e = mesh_cellgroup(:, 2);
    cellgroup_v = mesh_cellgroup(:, 3);
    cellgroup_dn = mesh_cellgroup(:, 4);
    cellgroup_de = mesh_cellgroup(:, 5);
    cellgroup_dv = mesh_cellgroup(:, 6);
end

modelsize2 = length(nonzeros(cellgroup_model));
disp(['Compression ratio of model cells is ' num2str(modelsize1/modelsize2) '.']);
 
end