%   #GravForModeling3D#
%   the algorithm GravForModeling3D provides a fast forward modeling of gravity data.
%   Copyright Tao Chen 2018/08/20
%   Please report any bug to geogoku@aliyun.com
%
%   Use this code, please refer to the paper 
%                       "Chen, T. and Zhang, G., 2018. Forward modeling of gravity anomalies
%                        based on cell mergence and parallel computing. Computers & Geosciences,
%                        120, 1-9. doi: 10.1016/j.cageo.2018.07.007."
%
%   The program is tested using Matlab 2014.

clc; 
close all; 
clear all; 
format compact;

%% Input parameters, please refer to UBC-MeshTools3d for the format of file_mesh and file_model.
meshFile = 'mesh32x32x32';             % file name of the mesh
modelFile = 'model32x32x32';          % file name of the model

inputFlag = 1;
% inputFlag:
%                    1 User specifies the input parameters of the survey grid.
%                    2 read parameters of the survey grid from the parameter file.
%                    3 program adaptive setting the parameters of the survey grid.

Mode = 4;
% Mode£¨par_flag ºÍ flag_merge£©£º
%                    1 Serial algorithm
%                    2 Parallel program (if the device supports GPU, par_flag = 3, else par_flag = 2)
%                    3 mesh merging + Serial algorithm
%                    4 mesh merging + Parallel program
% parFlag:
%                    1 serial program
%                    2 multicore parallel program
%                    3 GPU parallel program (auto detection)
% mergeFlag: 
%                    0 Do not execute mesh merging
%                    1 Execute mesh merging

%% argument setting
ngpus = gpuDeviceCount;              % determine if the system supports GPU operations

switch Mode
    case 1
        mergeFlag = 0;
        parFlag = 1;
    case 2
        mergeFlag = 0;
        if ngpus>0
            parFlag = 3;
        else
            parFlag = 2;
            MatlabParallelPool('mode', 'start');
        end
    case 3
        mergeFlag = 1;
        parFlag = 1;
    case 4
        mergeFlag = 1;
        if ngpus>0
            parFlag = 3;
        else
            parFlag = 2;
            MatlabParallelPool('mode', 'start');
        end
    otherwise
        error('GravFor3D: ProgramInput: Invalid_Mode');
end

%% Setup the preliminaries for the forward modeling
% survey grid info
tic;
switch inputFlag
    case 1                                                             % Specify input parameters of the survey grid
        RNx = 64;                                                 % Point number of receiver in x(north) direction
        RNy = 64;                                                 % Point number of receiver in y(east) direction
        RNz = 1;                                                     % Point number of receiver in z(up) direction
        Min_x = 0;                                                  % Minimum coordinate in x(north) direction
        Max_x = 100000;                                          % Maximun coordinate in x(north) direction
        Min_y = 0;                                                  % Minimum coordinate in y(east) direction
        Max_y = 100000;                                          % Maximun coordinate in y(east) direction
        Min_z = 1;                                                   % Minimum coordinate in z(up) direction
        Max_z = 1;                                                  % Maximun coordinate in z(up) direction
        if RNz == 1
            RSx = (Max_x - Min_x)/(RNx - 1);          % Point spacing of receiver in x(north) direction
            RSy = (Max_y - Min_y)/(RNy - 1);          % Point spacing of receiver in y(east) direction
            nd = RNy * RNx;                                    % Total number of the survey grid
            RPz = ones(nd, 1) * Min_z;
            i = 1 : RNx;      j = 1 : RNy;
            RPx = Min_x + (i - 1)' * RSx;
            RPy = Min_y + (j - 1)' * RSy;
            RPx = repmat(RPx, RNy, 1);                   % coordinates in x direction
            RPy = repmat(RPy', RNx, 1);                        
            RPy = reshape(RPy, nd, 1);                    % coordinates in y direction
            clear j RSy Max_z  Min_z Max_y Min_y;
        else
            RSx = (Max_x - Min_x)/(RNx - 1);          % Point spacing of receiver in x(north) direction
            RSy = (Max_y - Min_y)/(RNy - 1);          % Point spacing of receiver in y(east) direction
            RSz = (Max_z - Min_z)/(RNz - 1);          % Point spacing of receiver in z(up) direction
            nd = RNx * RNy * RNz;                         % Total number of the survey grid
            i = 1 : RNx; j = 1 : RNy; k = 1 : RNz;
            RPx = Min_x + (i - 1)' * RSx;
            RPy = Min_y + (j - 1)' * RSy;
            RPx = repmat(RPx, RNy, 1);                      
            RPy = repmat(RPy', RNx, 1);                        
            RPy = reshape(RPy, nd, 1);                       
            RPx = repmat(RPx, RNz, 1);                      % coordinates in x direction
            RPy = repmat(RPy, RNz, 1);                      % coordinates in y direction
            RPz = Min_z + (k - 1) * RSz;
            RPz = reshape(RPz, 1, 1, RNz);
            RPz = repmat(RPz, RNx, RNy, 1);
            RPz = reshape(RPz, nd, 1);                       % coordinates in z direction
        end
        clear j k ip RSy RSz Max_z Min_z Max_y Min_y Max_x Min_x RSx;
    case 2                                                % read parameters from the parameter file
        % read headerline info
        fp = fopen('obs.loc','r');
        headerline = fgets(fp);
        fclose(fp);
        NumExtracted = regexp(headerline,'\d*\.?\d*','match');
        nd = str2num(cell2mat(NumExtracted(1)));
        [RPy, RPx, RPz] = textread('obs.loc','%f%f%f','headerlines',1);
        clear fid headerline NumExtracted;
    case 3
        RPx = unique(meshcell(:, 1));
        RPy = unique(meshcell(:, 2));
        nn = length(RPx);  ne = length(RPy);
        RPx = repmat(RPx, ne, 1);                             % coordinates in x direction
        RPy = repmat(RPy', nn, 1);                        
        RPy = reshape(RPy, nn * ne, 1);                    % coordinates in y direction
        nd = nn * ne;
        RPz = RPy;
        if min(meshcell(:, 4)) > min(meshcell(:, 5))
            RPz(:) = min(meshcell(:, 5));
        else
            RPz(:) = min(meshcell(:, 4));
        end
end
time1 = toc;

%% mesh merge
tic;
if mergeFlag
    disp('Model mesh merging...');
    [cell_n, cell_e, cell_v, cell_dn, cell_de, cell_dv, model] = CellMerge3D(meshFile, modelFile, 0, 1);
else   % read model and mesh info
    [~, cell_n, cell_e, cell_v, cell_dn, cell_de, cell_dv] = Mesh3D(meshFile);
    model = importdata(modelFile);
end
time2 = toc;

%% caculate the anomaly data
disp('Computing of anomalies...');
tic;
anomaly_ds_gpu = GravFor3D(mergeFlag, parFlag, cell_n, cell_e, cell_v, cell_dn, cell_de, cell_dv, model, RPx, RPy, RPz);
time3 = toc;
clear mn nd mesh_info ngpus

%% Mapping
disp('Mapping...')

contourf(reshape(RPx, RNx, RNy), reshape(RPy, RNx, RNy), ... 
    reshape(anomaly_ds_gpu, RNx, RNy), 20); ...
    h = colorbar; set(get(h,'Title'),'string','mGal'); xlabel('East (m)'); ylabel('North (m)');

disp(['CPU time of the GravForModeling3D is ' num2str(time1 + time2 + time3) ' s.']);

if parFlag == 2
    MatlabParallelPool('mode', 'close');
end

clear i Mode par_flag currentFolder flag_merge file_model file_mesh flag_input