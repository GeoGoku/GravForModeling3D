function [] = MatlabParallelPool(varargin)
%MATLABPARALLELPOOL starts or closes the matlab parallel pool.
%   MATLABPARALLELPOOL performs 'start' or 'close' action of the
%   matlab parallel pool.
%
%   MatlabParallelPool('mode', 'start') starts the matlab
%   parallel pool using the default poolNum or currentPoolNum.
%   
%   MatlabParallelPool('mode', 'start', 'size', poolNum) starts the matlab
%   parallel pool, the parameter 'size' specifies the poolNum.
%
%   MatlabParallelPool('mode', 'close') closes the matlab parallel pool.
%
%   Example:
%     poolNum = 4£»
%     MatlabParallelPool('mode', 'start', 'size', poolNum) 
%
%   T. Chen 24-Sep-2016.
%   Copyright 2016 T. Chen.

% Check arguments
nbIn = nargin;
switch nbIn
    case {0, 1}
        error('MatlabParallelPool: FunctionInput: NotEnough_ArgNum');
    case {3}
        error('MatlabParallelPool: FunctionInput: Invalid_ArgNum');
    case {2, 4}
    otherwise
        error('MatlabParallelPool: FunctionInput: TooMany_ArgNum');
end

% Initialize the switchMode and poolNum
for k = 1 : 2 : nbIn-1
    switch varargin{k}
        case 'mode'
            switchMode = varargin{k + 1};
        case 'size'
            poolNum = varargin{k + 1};
    end
end

if strcmp(switchMode, 'start')   % start the matlab parallel pool
    % check if the parallel pool has been started.
    p = gcp('nocreate');
    if isempty(p)
        currentPoolNum = 0;
    else
        currentPoolNum = p.NumWorkers;
    end
    % configure the parallel pool with the argument user specified.
    if currentPoolNum == 0
        if nargin == 2
            parpool('local');
        else
            try
                parpool('local', poolNum);
            catch ce
                parpool;
                display(ce.message);
            end
        end
    else
        disp('parpool has started.');
        if nargin == 2
            disp(['parallel pool has been started with ', num2str(currentPoolNum), ' workers.']);
        else
            if currentPoolNum ~= poolNum
                MatlabParallelPool('mode', 'close') ;
                MatlabParallelPool('mode', 'start', 'size', poolNum);
            end
        end
    end
elseif strcmp(switchMode, 'close')   % close the matlab parallel pool
    poolobj = gcp('nocreate');
    delete(poolobj);
end

end