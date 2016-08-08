%% PDSCH Throughput Conformance Test for Single Antenna (TM1), Transmit Diversity (TM2), Open Loop (TM3) and Closed Loop (TM4/6) Spatial Multiplexing

% * TM1: Single antenna (Port 0)
% * TM2: Transmit diversity
% * TM3: Open loop codebook based precoding: Cyclic Delay Diversity (CDD)
% * TM4: Closed loop codebook based spatial multiplexing
% * TM6: Single layer closed loop codebook based spatial multiplexing
% 
%                                 Alloc    Mod               Coding  Number
% RMC   TM                         RBs    scheme  CellRefP    rate   layers 
% -------------------------------------------------------------------------
% R.0  Port0                        1 RB   16QAM   CellRefP=1  R=1/2    1
% R.1  Port0                        1 RB   16QAM   CellRefP=1  R=1/2    1
% R.2  Port0                       50 RB    QPSK   CellRefP=1  R=1/3    1
% R.3  Port0                       50 RB   16QAM   CellRefP=1  R=1/3    1
% R.10 TxDiversity|SpatialMux      50 RB    QPSK   CellRefP=2  R=1/3    2
% R.11 TxDiversity|SpatialMux|CDD  50 RB   16QAM   CellRefP=2  R=1/2    2
% R.12 TxDiversity                  6 RB    QPSK   CellRefP=4  R=1/3    4
% R.13 SpatialMux                  50 RB    QPSK   CellRefP=4  R=1/3    1
%

%% 
clear;
clc;
close all;

NFrames = 10;                % Number of frames
SNRIn = 20;%[-5.0:5.0:5.0];   % SNR range

txMode = 'TM6'; % TM1, TM2, TM3, TM4, TM6

enb = []; % clear enb
switch txMode
    case 'TM1' % based on R.3
        fprintf('\nTM1 - Single antenna (port 0)\n');
        enb.RC = 'R.3'; % Port 0, 50 RB, 16QAM, CellRefP=1, R=1/2        
        enb.NDLRB = 50;
        enb.CellRefP = 1;        
        enb.PDSCH.TxScheme = 'Port0';
        enb.PDSCH.Modulation = {'16QAM'};
    case 'TM2' % based on R.11
        fprintf('\nTM2 - Transmit diversity\n');
        enb.RC = 'R.11'; % TxDiversity, 50 RB, 16QAM, CellRefP=2, R=1/2
        enb.NDLRB = 50;
        enb.CellRefP = 2;
        enb.PDSCH.TxScheme = 'TxDiversity';
        enb.PDSCH.Modulation = {'16QAM'};
        enb.PDSCH.NLayers = 2;   
    case 'TM3' % based on R.11
        fprintf('\nTM3 - CDD\n');
        enb.RC = 'R.11'; % CDD, 50 RB, 16QAM, CellRefP=2, R=1/2        
        enb.NDLRB = 50;
        enb.CellRefP = 2;        
        enb.PDSCH.TxScheme = 'CDD';
        enb.PDSCH.Modulation = {'16QAM', '16QAM'};
        enb.PDSCH.NLayers = 2;
    case 'TM4' % based on R.11
        fprintf('\nTM4 - Codebook based spatial multiplexing\n');
        enb.RC = 'R.11'; % SpatialMux, 50 RB, 16QAM, CellRefP=2, R=1/2
        enb.NDLRB = 50;
        enb.CellRefP = 2;              
        enb.PDSCH.TxScheme = 'SpatialMux';
        enb.PDSCH.Modulation = {'16QAM', '16QAM'};
        enb.PDSCH.NLayers = 2;
        enb.PDSCH.CodebookSubset = '110000'; % No codebook restriction
        enb.PDSCH.Rho = -3;
    case 'TM6' % based on R.13
        fprintf('\nTM6 - Codebook based spatial multiplexing with single layer\n');
        enb.RC = 'R.13'; % SpatialMux, 50 RB, QPSK, CellRefP=4, R=1/3
        enb.NDLRB = 50;
        enb.CellRefP = 4;        
        enb.PDSCH.TxScheme = 'SpatialMux';
        enb.PDSCH.Modulation = {'QPSK'};
        enb.PDSCH.NLayers = 1;
        enb.PDSCH.CodebookSubset = ''; % No codebook restriction
    otherwise
        error('Non-supported TM')
end

% Set enb fields applying to all TMs
enb.TotSubframes = 1; % Generate one subframe at a time
enb.PDSCH.CSI = 'On'; % Soft bits are weighted by channel state information

%%
% Now we call lteRMCDL to populate all the fields in |enb| not specified
% above. This function makes use of the specified RC and the user provided
% fields. The |enb| output will include all the required fields not
% specified by the user, in addition to the specified fields.
enb = lteRMCDL(enb);

%%
% The output |enb| structure contains, amongst other fields, the transport
% block sizes and redundancy version sequence for each codeword subframe
% within a frame. These will be used later in the simulation.
rvSequence = enb.PDSCH.RVSeq;
trBlkSizes = enb.PDSCH.TrBlkSizes;

%%
% Extract the number of codewords |ncw|. This is the number of entries in
% the |enb.PDSCH.Modulation| field.

ncw = length(enb.PDSCH.Modulation);

%%
% Set the PMI delay for the closed-loop TMs (TM4 and TM6). This is the
% delay between a PMI being passed from UE to eNodeB as defined in
% TS36.101, Table 8.2.1.4.2-1 [ <#21 1> ]. 

pmiDelay = 8;

%%
% Next we print a summary of some of the more relevant simulation
% parameters. Check these values to make sure they are as expected. The
% code rate displayed can be useful to detect problems if manually
% specifying the transport block sizes. Typical values are 1/3, 1/2 and
% 3/4.

fprintf('\n-- Parameter summary: --------------------------------------------\n');
disp(['                 Transmission scheme: ' enb.PDSCH.TxScheme]);
disp(['  Number of downlink resource blocks: ' num2str(enb.NDLRB)]);
disp([' Number of allocated resource blocks: ' num2str(length(enb.PDSCH.PRBSet))]);
disp(['Cell-specific reference signal ports: ' num2str(enb.CellRefP)]);
disp(['                 Number of codewords: ' num2str(ncw)]);
disp(['                 Transmission layers: ' num2str(enb.PDSCH.NLayers)]);
disp(['               Modulation codeword 1: ' enb.PDSCH.Modulation{1}]);
disp(['    Transport block sizes codeword 1: ' num2str(enb.PDSCH.TrBlkSizes(1,:))]);

TBS = enb.PDSCH.TrBlkSizes(1,1);
dlschInfo = lteDLSCHInfo(TBS);
[~,indicesInfo] = ltePDSCHIndices(enb,enb.PDSCH,enb.PDSCH.PRBSet);
codeRate = double(dlschInfo.Bout)/double(indicesInfo.G(1));
disp(['                Code rate codeword 1: ' num2str(codeRate)]);
if size(enb.PDSCH.TrBlkSizes,1)==2
    disp(['               Modulation codeword 2: ' enb.PDSCH.Modulation{2}]);
    disp(['    Transport block sizes codeword 2: ' num2str(enb.PDSCH.TrBlkSizes(1,:))]);
    TBS = enb.PDSCH.TrBlkSizes(2,1);
    dlschInfo = lteDLSCHInfo(TBS);
    codeRate = double(dlschInfo.Bout)/double(indicesInfo.G(2));
    disp(['                Code rate codeword 2: ' num2str(codeRate)]);
end
disp('------------------------------------------------------------------');

%% Propagation Channel Model Configuration
% The structure |channel| contains the channel model configuration
% parameters.

channel.Seed = 10;                   % Random channel seed
channel.NRxAnts = 2;                 % 2 receive antennas
channel.DelayProfile = 'ETU';        % Delay profile
channel.DopplerFreq = 70;            % Doppler frequency
channel.MIMOCorrelation = 'Low';     % Multi-antenna correlation
channel.NTerms = 16;                 % Oscillators used in fading model
channel.ModelType = 'GMEDS';         % Rayleigh fading model type
channel.InitPhase = 'Random';        % Random initial phases
channel.NormalizePathGains = 'On';   % Normalize delay profile power  
channel.NormalizeTxAnts = 'On';      % Normalize for transmit antennas

% The sampling rate for the channel model is set using the value returned
% from <matlab:doc('lteOFDMInfo') lteOFDMInfo>.

ofdmInfo = lteOFDMInfo(enb);
channel.SamplingRate = ofdmInfo.SamplingRate;

%% Channel Estimator Configuration
% The variable |perfectChanEstimator| controls channel estimator behavior.
% Valid values are |true| or |false|. When set to |true| a perfect channel
% response is used as estimate, otherwise an imperfect estimation based on
% the values of received pilot signals is obtained.

% Perfect channel estimator flag
perfectChanEstimator = false;

%%
% If |perfectChanEstimator| is set to false a configuration structure |cec|
% is needed to parameterize the channel estimator. A Doppler frequency of
% 70 Hz causes the channel to fade relatively slowly over time, therefore a
% large time averaging window can be used. No averaging is performed in the
% frequency domain.

% Configure channel estimator
cec.PilotAverage = 'UserDefined';   % Type of pilot symbol averaging
cec.FreqWindow = 1;                 % Frequency window size in REs
cec.TimeWindow = 31;                % Time window size in REs
cec.InterpType = 'Cubic';           % 2D interpolation type
cec.InterpWindow = 'Centered';      % Interpolation window type
cec.InterpWinSize = 1;              % Interpolation window size


%% Setup HARQ Processes 
% This example models an implementation of a HARQ retransmission scheme.
% The HARQ table stores the index of each HARQ processes to use. This is
% generated for 8-HARQ processes using the helper function
% <matlab:edit('hHARQTable.m') hHARQTable.m>.

[harqTable, noHarqProcesses] = hHARQTable();   

%% Display Simulation Information
% The variable |displaySimulationInformation| controls the display of
% simulation information such as the HARQ process number used for each
% subframe. In case of CRC error the value of the index to the RV sequence
% is also displayed.

displaySimulationInformation = true;

%% Processing Loop

dims = lteDLResourceGridSize(enb);
P = dims(3);

% Initialize variables used in the simulation and analysis
% Array to store the maximum throughput for all SNR points
maxThroughput = zeros(length(SNRIn),1); 
% Array to store the simulation throughput for all SNR points
simThroughput = zeros(length(SNRIn),1);

% The temporary variables 'enb_init' and 'channel_init' are used to create
% the temporary variables 'enb' and 'channel' within the SNR loop to create
% independent simulation loops for the parfor loop
enb_init = enb;
channel_init = channel;
legendString = ['Throughput: ' enb.PDSCH.TxScheme];
allRvSeqPtrHistory = cell(1,numel(SNRIn));
nFFT = ofdmInfo.Nfft;

for snrIdx = 1:numel(SNRIn)

    rng(snrIdx,'combRecursive');
    
    SNRdB = SNRIn(snrIdx);
    fprintf('\nSimulating at %gdB SNR for %d Frame(s)\n' ,SNRdB, NFrames);
    
    % Initialize variables used in the simulation and analysis
    offsets = 0;            % Initialize frame offset value
    offset = 0;             % Initialize frame offset value for the radio frame
    blkCRC = [];            % Block CRC for all considered subframes
    bitTput = [];           % Number of successfully received bits per subframe
    txedTrBlkSizes = [];    % Number of transmitted bits per subframe
    enb = enb_init;         % Initialize RMC configuration
    channel = channel_init; % Initialize channel configuration
    pmiIdx = 0;             % PMI index in delay queue
    
    % The variable harqPtrTable stores the history of the value of the
    % pointer to the RV sequence values for all the HARQ processes
    rvSeqPtrHistory = zeros(ncw, NFrames*10);        
    
    % Initialize state of all HARQ processes
    harqprocess = hDownlinkNewHARQProcess(enb);
    harqProcesses = repmat(harqprocess, 1, max(harqTable));
    
    % Use random PMIs for the first 'pmiDelay' subframes until feedback is
    % available from the UE; note that PMI feedback is only applicable for
    % spatial multiplexing TMs (TM4 and TM6), but the code here is required
    % for complete initialization of variables in the SNR loop when using
    % the Parallel Computing Toolbox.
    pmidims = ltePMIInfo(enb,enb.PDSCH);
    txPMIs = randi([0 pmidims.MaxPMI], pmidims.NSubbands, pmiDelay);

    for subframeNo = 0:(NFrames*10-1)
        
        % Reinitialize channel seed for each subframe to increase
        % variability
        channel.Seed = 1+subframeNo;
        
        % Update subframe number
        enb.NSubframe = subframeNo;

        % Get HARQ index for given subframe from HARQ index table
        harqIdx = harqTable(mod(subframeNo, length(harqTable))+1);

        % Update current HARQ process
        harqProcesses(harqIdx) = hHARQScheduling( ...
            harqProcesses(harqIdx), subframeNo, rvSequence);
        
        % Extract the current subframe transport block size(s)
        trBlk = trBlkSizes(:, mod(subframeNo, 10)+1).';
        
        % Display run time information
        if displaySimulationInformation && (any(trBlk) ~= 0)
            % map harqIdx to harq process number
            harqProcNo = find([1:5 7:9]==harqIdx);
            disp(' ');
            disp(['Subframe: ' num2str(subframeNo) '. HARQ process index: ' num2str(harqProcNo)]);
        end
        
        % Update RV sequence pointer table
        if any(trBlk) ~= 0
            rvSeqPtrHistory(:,subframeNo+1) = harqProcesses(harqIdx).txConfig.RVIdx.';
        else
            rvSeqPtrHistory(:,subframeNo+1) = NaN; % no data in subframe
        end

        % Update the PDSCH transmission config with HARQ process state
        enb.PDSCH = harqProcesses(harqIdx).txConfig;      
        data = harqProcesses(harqIdx).data;

        % Set the PMI to the appropriate value in the delay queue
        if strcmpi(enb.PDSCH.TxScheme,'SpatialMux')
            pmiIdx = mod(subframeNo, pmiDelay);     % PMI index in delay queue
            enb.PDSCH.PMISet = txPMIs(:, pmiIdx+1); % Set PMI
        end

        % Create transmit waveform and add 25 sample padding. This is to
        % cover the range of delays expected from channel modeling (a
        % combination of implementation delay and channel delay spread)
        txWaveform = [lteRMCDLTool(enb, data); zeros(25, P)];

        % Initialize channel time for each subframe
        channel.InitTime = subframeNo/1000;

        % Pass data through channel model
        rxWaveform = lteFadingChannel(channel, txWaveform);

        % Calculate noise gain including compensation for downlink power
        % allocation
        SNR = 10^((SNRdB-enb.PDSCH.Rho)/20);

        % Normalize noise power to take account of sampling rate, which is
        % a function of the IFFT size used in OFDM modulation, and the 
        % number of antennas
        N0 = 1/(sqrt(2.0*enb.CellRefP*double(nFFT))*SNR);

        % Create additive white Gaussian noise
        noise = N0*complex(randn(size(rxWaveform)), ...
                            randn(size(rxWaveform)));

        % Add AWGN to the received time domain waveform        
        rxWaveform = rxWaveform + noise;

        % Once every frame, on subframe 0, calculate a new synchronization
        % offset
        if (mod(subframeNo,10) == 0)
            offset = lteDLFrameOffset(enb, rxWaveform);
            if (offset > 25)
                offset = offsets(end);
            end
            offsets = [offsets offset]; %#ok
        end
        
        % Synchronize the received waveform
        rxWaveform = rxWaveform(1+offset:end, :);

        % Perform OFDM demodulation on the received data to recreate the
        % resource grid
        rxSubframe = lteOFDMDemodulate(enb, rxWaveform);

        % Channel estimation
        if(perfectChanEstimator) 
            estChannelGrid = lteDLPerfectChannelEstimate(enb, channel, offset); %#ok
            noiseGrid = lteOFDMDemodulate(enb, noise(1+offset:end ,:));
            noiseEst = var(noiseGrid(:));
        else
            [estChannelGrid, noiseEst] = lteDLChannelEstimate( ...
                enb, enb.PDSCH, cec, rxSubframe);
        end

        % Get PDSCH indices
        [pdschIndices,~] = ...
            ltePDSCHIndices(enb, enb.PDSCH, enb.PDSCH.PRBSet);
        
        % Get PDSCH resource elements from the received subframe. Scale the
        % received subframe by the PDSCH power factor Rho. The PDSCH is
        % scaled by this amount, while the cell reference symbols used for
        % channel estimation (used in the PDSCH decoding stage) are not.
        [pdschRx, pdschHest] = lteExtractResources(pdschIndices, ...
            rxSubframe*(10^(-enb.PDSCH.Rho/20)), estChannelGrid);
        
        % Decode PDSCH
        [dlschBits,~] = ltePDSCHDecode(enb, enb.PDSCH,...
            pdschRx, pdschHest, noiseEst);
        
        % Decode the DL-SCH
        [decbits, harqProcesses(harqIdx).blkerr, ...
            harqProcesses(harqIdx).decState] = lteDLSCHDecode( ...
            enb, enb.PDSCH, trBlk, dlschBits, ...
            harqProcesses(harqIdx).decState);
        
        % Display block errors
        if displaySimulationInformation && any(harqProcesses(harqIdx).blkerr)
            disp(['Block error. RV index: ' num2str(harqProcesses(harqIdx).txConfig.RVIdx) ', CRC: ' num2str(harqProcesses(harqIdx).blkerr)])
        end
        
        % Store values to calculate throughput
        % Only for subframes with data
        if(any(trBlk) ~= 0)
            blkCRC = [blkCRC harqProcesses(harqIdx).blkerr]; %#ok<AGROW>
            bitTput = [bitTput trBlk.*(1- ...
                harqProcesses(harqIdx).blkerr)]; %#ok<AGROW>
            txedTrBlkSizes = [txedTrBlkSizes trBlk]; %#ok<AGROW>
        end

        % Provide PMI feedback to the eNodeB
        if strcmpi(enb.PDSCH.TxScheme,'SpatialMux')
            PMI = ltePMISelect(enb, enb.PDSCH, estChannelGrid, noiseEst);
            txPMIs(:, pmiIdx+1) = PMI;
        end
    end
    
    % Calculate maximum and simulated throughput
    maxThroughput(snrIdx) = sum(txedTrBlkSizes); % Maximum possible throughput
    simThroughput(snrIdx) = sum(bitTput,2);      % Simulated throughput
    
    % Display the results dynamically in the command window
    fprintf([['\nThroughput(Mbps) for ', num2str(NFrames) ' Frame(s) '],...
        '= %.4f\n'], 1e-6*simThroughput(snrIdx)/(NFrames*10e-3));
    fprintf(['Throughput(%%) for ', num2str(NFrames) ' Frame(s) = %.4f\n'],...
        simThroughput(snrIdx)*100/maxThroughput(snrIdx));
    
    allRvSeqPtrHistory{snrIdx} = rvSeqPtrHistory;
    
end

% Plot the RV sequence for all HARQ processes
for snrIdx = 1:numel(SNRIn)
    SNRdB = SNRIn(snrIdx);
    rvSeqPtrHistory = allRvSeqPtrHistory{snrIdx};
    figure;
    line((0:(NFrames*10-1)),rvSeqPtrHistory(1,:),'Color','k',...
        'Marker','o')
    ax1 = gca; % current axes
    ax1.YLabel.String = 'cw1  RV sequence index';
    ax1.XLabel.String = 'subframe';
    ax1.YLim = [min(rvSeqPtrHistory(:))-1 max(rvSeqPtrHistory(:))+1];
    if ncw == 2
        ax1_pos = ax1.Position; % position of first axes
        ax2 = axes('Position',ax1_pos,...
            'YAxisLocation','right','Color','none','YColor','r');
        ax2.YLabel.String = 'cw2  RV sequence index';
        ax2.YLim = ax1.YLim;
        line((0:(NFrames*10-1)),rvSeqPtrHistory(2,:),'Parent',ax2,...
            'Color','r','Marker','+','LineStyle',':')
    end
    title(['RV Sequence Pointers: SNR ' num2str(SNRdB) ' dB'])
end

%% RV Sequence Pointer Plots
% The code above also generates plots with the value of the pointers to the
% elements in the RV sequence for the simulated subframes. This provides an
% idea of the retransmissions required. We plot the pointers and note the
% RV values used in case these are not organized in ascending order. For
% example, in some cases the RV sequence can be [0, 2, 3, 1]. Plotting
% these values as they are used will not provide a clear idea of the number
% of retransmissions needed.
% 
% When transmitting a new packet, the first element of the RV sequence is
% used. In the plots above a value of 1 is shown for that subframe. This is
% the case at the beginning of the simulation. If a retransmission is
% required, the next element in the RV sequence is selected and the pointer
% is increased. A value of 2 will be plotted for the subframe where the
% retransmission takes place. If further retransmissions are required the
% pointer value will increase further. Note that the plots do not show any
% value in subframe 5 of consecutive frames. This is because no data is
% transmitted in those subframes.
%
% The figure shown below was obtained simulating 10 frames. Note how in
% some cases up to 3 retransmissions are required.
%
% <<10FramesRVSeqPointer.jpg>>

%% Throughput Results
% The throughput results for the simulation are displayed in the MATLAB(R)
% command window after each SNR point is completed. They are also captured
% in |simThroughput| and |maxThroughput|. |simThroughput| is an array with
% the measured throughput in number of bits for all simulated SNR points.
% |maxThroughput| stores the maximum possible throughput in number of bits
% for each simulated SNR point.

% Plot throughput
figure
plot(SNRIn, simThroughput*100./maxThroughput,'*-.');
xlabel('SNR (dB)');
ylabel('Throughput (%)');
legend(legendString,'Location','NorthWest');
grid on;

%%
% The generated plot has been obtained with a low number of frames,
% therefore the results shown are not representative. A longer simulation
% obtained with 1000 frames produced the results shown below.
%
% <<TM4throughout1000frames.jpg>>

%% Appendix
% This example uses the helper functions:
%
% * <matlab:edit('hHARQTable.m') hHARQTable.m>
% * <matlab:edit('hHARQScheduling.m') hHARQScheduling.m>
% * <matlab:edit('hDownlinkNewHARQProcess.m') hDownlinkNewHARQProcess.m>

%% Selected Bibliography
% # 3GPP TS 36.101 "User Equipment (UE) radio transmission and reception"

displayEndOfDemoMessage(mfilename)
