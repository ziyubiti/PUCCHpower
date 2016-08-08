%% PUCCH2 CQI BLER 

clc;
clear;
close all;
path(path,'../Model');

numSubframes = 10000;                          % Number of subframes
SNRdB = [-15:1.0:-2.0];         % SNR range
NTxAnts = 1;                                % Number of transmit antennas

%% UE Configuration
ue = struct;                                % UE config structure
ue.NULRB = 100;                               % 6 resource blocks
ue.CyclicPrefixUL = 'Normal';               % Normal cyclic prefix
ue.Hopping = 'Off';                         % No frequency hopping
ue.NCellID = 9;
ue.RNTI = 1;                                % Radio network temporary id
ue.NTxAnts = NTxAnts;

%%   PUCCH 2 Configuration                  
ACK = [];
pucch = struct; % PUCCH config structure
% Vector of PUCCH resource indices, one per transmission antenna. This is
% the n2pucch parameter
pucch.ResourceIdx = 0:ue.NTxAnts-1;
% Set the size of resources allocated to PUCCH format 2
pucch.ResourceSize = 0;
% Number of cyclic shifts used for PUCCH format 1 in resource blocks with a
% mixture of formats 1 and 2. This is the N1cs parameter
pucch.CyclicShifts = 0;   

%% Propagation Channel Configuration
channel = struct;                   % Channel config structure
channel.NRxAnts = 2;                % Number of receive antennas
channel.DelayProfile = 'EVA';       % Channel delay profile
channel.DopplerFreq = 5.0;         % Doppler frequency in Hz
channel.MIMOCorrelation = 'Low';    % Low MIMO correlation
channel.NTerms = 16;                % Oscillators used in fading model
channel.ModelType = 'GMEDS';        % Rayleigh fading model type    
channel.Seed = 3;                   % Random number generator seed 
channel.InitPhase = 'Random';       % Random initial phases     
channel.NormalizePathGains = 'On';  % Normalize delay profile power   
channel.NormalizeTxAnts = 'On';     % Normalize for transmit antennas

% SC-FDMA modulation information: required to get the sampling rate
info = lteSCFDMAInfo(ue);
channel.SamplingRate = info.SamplingRate;   % Channel sampling rate

%%  Channel Estimator Configuration

cec = struct;                     % Channel estimation config structure
cec.PilotAverage = 'UserDefined'; % Type of pilot averaging
cec.FreqWindow = 12;              % Frequency averaging window in REs (special mode)
cec.TimeWindow = 1;               % Time averaging window in REs (Special mode)     
cec.InterpType = 'cubic';         % Cubic interpolation

%% Simulation Loop for Configured SNR Points

BLER = zeros(size(SNRdB));
for nSNR = 1:length(SNRdB)

    % Detection failures counter
    failCount = 0;
 
    SNR = 10^(SNRdB(nSNR)/20);              % Convert dB to linear 
    N = 1/(SNR*sqrt(double(info.Nfft)))/sqrt(2.0*ue.NTxAnts); 
    rng('default');
    fprintf('\nSimulating at %g dB SNR for a total %d SubFrame(s)', ...
        SNRdB(nSNR), numSubframes);
    
    % Loop for subframes    
    offsetused = 0;    
    for nsf = 1:numSubframes        

        % Create resource grid        
        ue.NSubframe = mod(nsf-1, 10);   % Subframe number       
        reGrid = lteULResourceGrid(ue);  % Resource grid

        % Create PUCCH 2 and its DRS
        CQI = randi([0 1], 4, 1);             % Generate 4 CQI bits to send
        % Encode CQI bits to produce 20 bits
        coded = lteUCIEncode(CQI);              
        pucch2Sym = ltePUCCH2(ue, pucch, coded);     % PUCCH 2 modulation
        pucch2DRSSym = ltePUCCH2DRS(ue, pucch, ACK); % PUCCH 2 DRS creation

        % Generate indices for PUCCH 2 and its DRS      
        pucch2Indices = ltePUCCH2Indices(ue, pucch);
        pucch2DRSIndices = ltePUCCH2DRSIndices(ue, pucch);

        % Map PUCCH 2 and its DRS to the resource grid
        reGrid(pucch2Indices) = pucch2Sym;
        reGrid(pucch2DRSIndices) = pucch2DRSSym;

        % SC-FDMA modulation
        txwave = lteSCFDMAModulate(ue, reGrid);  
        channel.InitTime = (nsf-1)/1000;
        rxwave = lteFadingChannel(channel, [txwave;zeros(25, ue.NTxAnts)]);

        % Add noise at receiver
        noise = N*complex(randn(size(rxwave)), randn(size(rxwave)));
        rxwave = rxwave + noise;

        % Receiver 

        % Synchronization
        [offset, rxACK] = lteULFrameOffsetPUCCH2( ...
            ue, pucch, rxwave, length(ACK));
        if (offset<25)
            offsetused = offset;
        end

        % SC-FDMA demodulation
        rxgrid = lteSCFDMADemodulate(ue, rxwave(1+offsetused:end, :));

        % Channel estimation            
        [H, n0] = lteULChannelEstimatePUCCH2(ue, pucch, cec, rxgrid, rxACK);

        % Extract REs corresponding to the PUCCH 2 from the given subframe
        % across all receive antennas and channel estimates
        [pucch2Rx, pucch2H] = lteExtractResources(pucch2Indices, rxgrid, H);

        % MMSE Equalization
        eqgrid = lteULResourceGrid(ue);    
        eqgrid(pucch2Indices) = lteEqualizeMMSE(pucch2Rx, pucch2H, n0);      
        
        % PUCCH 2 demodulation
        rxBits = ltePUCCH2Decode(ue, pucch, eqgrid(pucch2Indices));

        % PUCCH 2 decoding
        decoded = lteUCIDecode(rxBits, length(CQI));                       

        % Record any decoding failures
        if (sum(decoded~=CQI)~=0)          
            failCount = failCount + 1;
        end

        % Perform PUCCH 2 DRS decoding. This is not required as part of
        % this test, but illustrates the steps involved.

        % Extract REs corresponding to the PUCCH 2 DRS from the given
        % subframe across all receive antennas and channel estimates
        [drsRx, drsH] = lteExtractResources(pucch2DRSIndices, rxgrid, H);

        % PUCCH 2 DRS Equalization
        eqgrid(pucch2DRSIndices) = lteEqualizeMMSE(drsRx, drsH, n0); 

        % PUCCH 2 DRS decoding
        rxACK = ltePUCCH2DRSDecode( ...
            ue, pucch, length(ACK), eqgrid(pucch2DRSIndices));

    end
    
    % Probability of erroneous block detection
    BLER(nSNR) = (failCount/numSubframes);

end

%% Results

semilogy(SNRdB, BLER, 'b-o', 'LineWidth', 2, 'MarkerSize', 7); 
grid on;    
% plot(-4.4, 0.01, 'rx', 'LineWidth', 2, 'MarkerSize', 7);
xlabel('SNR (dB)');
ylabel('CQI BLER');
title(['Wide CQI detection for ', num2str(numSubframes), ' SubFrames']);   
axis([SNRdB(1)-0.1 SNRdB(end)+0.1 0.0001 1.0]);
legend('PUCCH format 2');

displayEndOfDemoMessage(mfilename)
