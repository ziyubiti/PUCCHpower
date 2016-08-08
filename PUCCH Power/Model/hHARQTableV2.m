%hHARQTable HARQ process schedule
%   HARQTABLE = hHARQTable() is the HARQ process schedule for 8 processes.
%
%   HARQTABLE is an array of indices which are used to index a HARQ process
%   on a subframe-by-subframe basis. Therefore HARQTABLE provides process
%   indices for an integer number of frames. If a HARQ process cannot be
%   used, -1 is used.

%   Copyright 2009-2014 The MathWorks, Inc.

function harqTable = hHARQTableV2(DupMode)
    if strcmp(DupMode,'FDD') == 1
        noHarqProcesses = 8;  % Set number of HARQ processes, 8 for FDD,   7 for TDD FRC
    else
        noHarqProcesses = 7;
    end;
    
    % Create empty schedule over 10 subframes and as many frames
    % required until the schedule repeats
    harqTable = zeros(10,lcm(8,noHarqProcesses-1)/8);

    % Populate schedule
    harqTable(1, :) = 1;  % Subframe #0 always uses process #1
    harqTable(6, :) = -1; % Subframe #5 carries no data
    restTable = repmat(2:noHarqProcesses, ...
        1, lcm(8, noHarqProcesses-1)/(noHarqProcesses-1));
    harqTable(harqTable==0) = restTable;

    % As subframe #5 carries no data but still requires a HARQ process
    % index treat process #6 as unused and increment remaining
    % processes.
    ind = harqTable>5;
    harqTable(ind) = mod(harqTable(ind)+1, 10);
    harqTable(harqTable==-1) = 6;
    
    % Linearize table for subframe-by-subframe indexing
    harqTable = harqTable(:).';

end