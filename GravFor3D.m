function [gravAnomaly] = GravFor3D(mergeFlag, parFlag, xm, ym, zm, dx, dy, dz, model, RPx, RPy, RPz)
%GRAVFOR3D calculates the gravity anomaly using input parameters.
%   GRAVFOR3D(mergeFlag, parFlag, xm, ym, zm, dx, dy, dz, model, RPx, RPy, RPz)
%   calculates the gravity anomalies of the model on the survey grid.
%   The model info is contained in 'xm', 'ym', 'zm', 'dx', 'dy', 'dz' and 'model',
%   The info of the survey grid is contained in 'RPx', 'RPy' and 'RPz'.
%   User can change the parameters 'mergeFlag' and 'parFlag' to speed up
%   the program.
%
%   Input£º    
%          mergeFlag           flag value of merge mode
%          parFlag                flag value of parallel computing mode(1 for serial, 2 for CPU parallel, 3 for GPU parallel)
%          xm                       center coordinates of mesh cells in x direction (north is positive)
%          ym                       center coordinates of mesh cells in y direction (east is positive)
%          zm                       center coordinates of mesh cells in z direction (vertical down is positive)
%          dx                        size of mesh cells in x direction
%          dy                        size of mesh cells in y direction
%          dz                        size of mesh cells in z direction
%          model                  physical parameters of the model
%          RPx                      coordinates of survey grid in x direction 
%          RPy                      coordinates of survey grid in y direction 
%          RPz                      coordinates of survey grid in z direction 
%
%   Output£º
%          gravAnomaly       gravity anomaly
%
%   T. Chen 20-May-2016.
%   Copyright 2016-2018 T. Chen.

nbIn = nargin;
switch nbIn
    case {12}
    case {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}
        error('GravFor3D: FunctionInput: NotEnough_ArgNum');
    otherwise
        error('GravFor3D: FunctionInput: TooMany_ArgNum');
end

G0 = 6.67428e-11;   % gravitational constant
nd = length(RPy);    % number of survey data

 % mesh info
if mergeFlag == 0
    index = find(model);
    xm = xm(index);
    ym = ym(index);
    zm = zm(index);
    model = model(index);
    dx = dx(index);
    dy = dy(index);
    dz = dz(index);
    mn = length(index);   % number of mesh cells
else
    mn = length(model);   % number of mesh cells
end

switch parFlag
    case 1   % serial program
        gravAnomaly = zeros(nd, 1);
        if mn > nd
            for ip = 1 : nd
                x1 = xm - dx/2.0 - RPx(ip);
                x2 = x1 + dx;
                y1 = ym - dy/2.0 - RPy(ip);
                y2 = y1 + dy;
                z1 = zm - dz/2.0 + RPz(ip);
                z2 = z1 + dz;
                [kenelValue] = arrayfun(@GravKernel3D, x1, x2, y1, y2, z1, z2);  
                G = - G0 * 1.0e+8 * kenelValue;
                gravAnomaly(ip) = G' * model;
            end
        else
            for ct = 1 : mn
                x1 = xm(ct) - dx(ct)/2.0 - RPx;
                x2 = x1 + dx(ct);
                y1 = ym(ct) - dy(ct)/2.0 - RPy;
                y2 = y1 + dy(ct);
                z1 = zm(ct) - dz(ct)/2.0 + RPz;
                z2 = z1 + dz(ct);
                [kenelValue] = arrayfun(@GravKernel3D, x1, x2, y1, y2, z1, z2);  
                G = - G0 * 1.0e+8 * kenelValue;
                pointAnomaly = G * model(ct);
                gravAnomaly = gravAnomaly + pointAnomaly;
            end
        end
    case 2                          % multi CPU cores parallel computing
        gravAnomaly = zeros(nd, 1);
        if mn > nd
            parfor ip = 1 : nd
                x1 = xm - dx/2.0 - RPx(ip);
                x2 = x1 + dx;
                y1 = ym - dy/2.0 - RPy(ip);
                y2 = y1 + dy;
                z1 = zm - dz/2.0 + RPz(ip);
                z2 = z1 + dz;
                [kenelValue] = arrayfun(@GravKernel3D, x1, x2, y1, y2, z1, z2);  
                G = - G0 * 1.0e+8 * kenelValue;
                gravAnomaly(ip) = G' * model;
            end
        else
            parfor ct = 1 : mn
                x1 = xm(ct) - dx(ct)/2.0 - RPx;
                x2 = x1 + dx(ct);
                y1 = ym(ct) - dy(ct)/2.0 - RPy;
                y2 = y1 + dy(ct);
                z1 = zm(ct) - dz(ct)/2.0 + RPz;
                z2 = z1 + dz(ct);
                [kenelValue] = arrayfun(@GravKernel3D, x1, x2, y1, y2, z1, z2);  
                G = - G0 * 1.0e+8 * kenelValue;
                pointAnomaly = G * model(ct);
                gravAnomaly = gravAnomaly + pointAnomaly;
            end
        end
    case 3                          % GPU parallel computing
        gravAnomaly = zeros(nd, 1, 'gpuArray');
        xm =  gpuArray(xm);
        dx =  gpuArray(dx);
        RPx = gpuArray(RPx);
        ym =  gpuArray(ym);
        dy =  gpuArray(dy);
        RPy = gpuArray(RPy);
        zm = gpuArray(zm);
        dz = gpuArray(dz);
        RPz = gpuArray(RPz);
        model = gpuArray(model);
        if mn > nd
            for ip = 1 : nd
                x1 = xm - dx/2.0 - RPx(ip);
                x2 = x1 + dx;
                y1 = ym - dy/2.0 - RPy(ip);
                y2 = y1 + dy;
                z1 = zm - dz/2.0 + RPz(ip);
                z2 = z1 + dz;
                [kenelValue] = arrayfun(@GravKernel3D, x1, x2, y1, y2, z1, z2); 
                G = - G0 * 1.0e+8 * kenelValue;
                gravAnomaly(ip) = G' * model;
            end
        else
            for ct = 1 : mn
                x1 = xm(ct) - dx(ct)/2.0 - RPx;
                x2 = x1 + dx(ct);
                y1 = ym(ct) - dy(ct)/2.0 - RPy;
                y2 = y1 + dy(ct);
                z1 = zm(ct) - dz(ct)/2.0 + RPz;
                z2 = z1 + dz(ct);
                [kenelValue] = arrayfun(@GravKernel3D, x1, x2, y1, y2, z1, z2);  
                G = - G0 * 1.0e+8 * kenelValue;
                pointAnomaly = G * model(ct);
                gravAnomaly = gravAnomaly + pointAnomaly;
            end
        end
end

end

function [gravEffect] = GravKernel3D(x1, x2, y1, y2, z1, z2)
% GRAVKERNEL3D calculates the gravity effect of a cell on a survey point
% Input parameters contain relative distance between the cell and the survey point
%
%   Reference: Li, X., Chouteau, M., 1998. Three-dimensional gravity modeling in all space. Surv.
%                     Geophys. 19 (4), 339¨C368.
%
%   Input:
%          relative distance between the cell and the survey point, 
%          there are: x1, x2, y1, y2, z1, z2
%
%   Output:
%          gravEffect          gravity effect
%
%   T. Chen 20-May-2016.
%   Copyright 2016-2018 T. Chen.

r1 = sqrt(x1 .* x1 + y1 .* y1 + z1 .* z1);
r2 = sqrt(x1 .* x1 + y1 .* y1 + z2 .* z2);
r3 = sqrt(x1 .* x1 + y2 .* y2 + z1 .* z1);
r4 = sqrt(x1 .* x1 + y2 .* y2 + z2 .* z2);
r5 = sqrt(x2 .* x2 + y1 .* y1 + z1 .* z1);
r6 = sqrt(x2 .* x2 + y1 .* y1 + z2 .* z2);
r7 = sqrt(x2 .* x2 + y2 .* y2 + z1 .* z1);
r8 = sqrt(x2 .* x2 + y2 .* y2 + z2 .* z2);

% equation (8) in Li (1998) 
gravEffect = - (x1 .* log(y1 + r1) + y1 .* log(x1 + r1) - z1 .* atan2(x1 .* y1, z1 .* r1))...
                      + (x1 .* log(y1 + r2) + y1 .* log(x1 + r2) - z2 .* atan2(x1 .* y1, z2 .* r2))...
                      + (x1 .* log(y2 + r3) + y2 .* log(x1 + r3) - z1 .* atan2(x1 .* y2, z1 .* r3))...
                       - (x1 .* log(y2 + r4) + y2 .* log(x1 + r4) - z2 .* atan2(x1 .* y2, z2 .* r4))...
                      + (x2 .* log(y1 + r5) + y1 .* log(x2 + r5) - z1 .* atan2(x2 .* y1, z1 .* r5))...
                       - (x2 .* log(y1 + r6) + y1 .* log(x2 + r6) - z2 .* atan2(x2 .* y1, z2 .* r6))...
                       - (x2 .* log(y2 + r7) + y2 .* log(x2 + r7) - z1 .* atan2(x2 .* y2, z1 .* r7))...
                      + (x2 .* log(y2 + r8) + y2 .* log(x2 + r8) - z2 .* atan2(x2 .* y2, z2 .* r8));      

end