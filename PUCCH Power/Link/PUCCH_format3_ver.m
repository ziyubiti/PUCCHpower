%% PUCCH3 ACK Missed Detection Probability Conformance Test
% This example measures the ACK missed detection probability using the
% LTE System Toolbox(TM) under the single user Physical Uplink Control
% Channel (PUCCH3) conformance test conditions as defined in TS 36.104
% Section 8.3.6.1 [ <#9 1> ].

% Copyright 2011-2014 The MathWorks, Inc.

%% Introduction
% This example uses a simulation length of 10 subframes. This value has
% been chosen to speed up the simulation. A larger value should be chosen
% to obtain more accurate results. The target defined in TS36.104 Section
% 8.3.6.1 [ <#9 1> ] for 10 MHz bandwidth (50 resource blocks) and a
% single transmit antenna is an Acknowledgment (ACK) missed detection
% probability not exceeding 1% at an SNR of -3.7 dB. The test is defined
% for 1 transmit antenna.
numSubframes = 10;                          % Number of subframes
SNRdB = [-9.7 -7.7 -5.7 -3.7 -1.7];         % SNR range
NTxAnts = 1;                                % Number of transmit antennas

%% UE configuration

ue = struct;                                % UE config structure
ue.NULRB = 50;                              % 50 resource blocks (10 MHz)
ue.CyclicPrefixUL = 'Normal';               % Normal cyclic prefix
ue.NTxAnts = NTxAnts;
ue.NCellID = 9;
ue.RNTI = 1;                                % Radio network temporary id
ue.Hopping = 'Off';                         % No frequency hopping
ue.Shortened = 0;                           % No SRS transmission

%% PUCCH 3 Configuration                    

% Vector of PUCCH resource indices, one per transmission antenna. This is
% the n3pucch parameter.
pucch = struct;
pucch.ResourceIdx = 0:ue.NTxAnts-1;

%% Propagation Channel Configuration
% Configure the channel model with the parameters specified in the tests
% described in TS36.104 Section 8.3.6.1 [ <#9 1> ].

channel = struct;                 % Channel config structure
channel.NRxAnts = 2;              % Number of receive antennas
channel.DelayProfile = 'EPA';     % Channel delay profile
channel.DopplerFreq = 5.0;        % Doppler frequency in Hz
channel.MIMOCorrelation = 'Low';  % Low MIMO correlation
channel.NTerms = 16;              % Oscillators used in fading model
channel.ModelType = 'GMEDS';      % Rayleigh fading model type    
channel.Seed = 4;                 % Random number generator seed  
channel.InitPhase = 'Random';     % Random initial phases         
channel.NormalizePathGains = 'On';% Normalize delay profile power  
channel.NormalizeTxAnts = 'On';   % Normalize for transmit antennas

% SC-FDMA modulation information: required to get the sampling rate
info = lteSCFDMAInfo(ue);
channel.SamplingRate = info.SamplingRate; 

%%  Channel Estimator Configuration
% The channel estimator is configured using a structure |cec|. Here cubic
% interpolation will be used with an averaging window of 12-by-1 Resource
% Elements (REs). This configures the channel estimator to use a special
% mode which ensures the ability to despread and orthogonalize the
% different overlapping PUCCH transmissions.

cec = struct;                     % Channel estimation config structure
cec.PilotAverage = 'UserDefined'; % Type of pilot averaging
cec.FreqWindow = 12;              % Frequency averaging window in REs (special mode)
cec.TimeWindow = 1;               % Time averaging window in REs (Special mode)     
cec.InterpType = 'cubic';         % Cubic interpolation

%% Simulation Loop for Configured SNR Points
% For each SNR point the loop below calculates the probability of
% successful ACK detection using information obtained from |NSubframes|
% consecutive subframes. The following operations are performed for each
% subframe and SNR values:
%
% * Create an empty resource grid
% * Generate and map PUCCH 3 and its Demodulation Reference Signal (DRS) to
% the resource grid
% * Apply SC-FDMA modulation
% * Send the modulated signal through the channel 
% * Receiver synchronization
% * Apply SC-FDMA demodulation
% * Estimate the channel
% * Minimum Mean Squared Error (MMSE) equalization
% * PUCCH 3 demodulation/decoding
% * Record decoding failures

% Preallocate memory for missed detection probability vector
PMISS = zeros(size(SNRdB));
for nSNR = 1:length(SNRdB)
    
    % Detection failures counter
    missCount = 0;
    falseCount = 0;

    % Noise configuration
    SNR = 10^(SNRdB(nSNR)/20);              % Convert dB to linear
    % The noise added before SC-FDMA demodulation will be amplified by the
    % IFFT. The amplification is the square root of the size of the IFFT.
    % To achieve the desired SNR after demodulation the noise power is
    % normalized by this value. In addition, because real and imaginary
    % parts of the noise are created separately before being combined into
    % complex additive white Gaussian noise, the noise amplitude must be
    % scaled by 1/sqrt(2*ue.NTxAnts) so the generated noise power is 1.
    N = 1/(SNR*sqrt(double(info.Nfft)))/sqrt(2.0*ue.NTxAnts);
    % Set the type of random number generator and its seed to the default
    % value
    rng('default'); 

    % Loop for subframes    
    offsetused = 0;
    for nsf = 1:numSubframes   

        % Create resource grid        
        ue.NSubframe = mod(nsf-1,10);       % Subframe number
        reGrid = lteULResourceGrid(ue);       % Resource grid

        % Generate PUCCH 3 and its DRS
        N_ACK = 4;                          % 4 ACK bits
        ACK = randi([0 1], N_ACK, 1);       % Generate N_ACK random bits
        coded = lteUCI3Encode(ACK);         % PUCCH 3 coding
        pucch3Sym = ltePUCCH3(ue, pucch, coded);     % PUCCH 3 modulation  
        pucch3DRSSym = ltePUCCH3DRS(ue, pucch, ACK); % PUCCH 3 DRS creation

        % Generate indices for PUCCH 3 and its DRS      
        pucch3Indices = ltePUCCH3Indices(ue, pucch);
        pucch3DRSIndices = ltePUCCH3DRSIndices(ue, pucch);

        % Map PUCCH 3 and its DRS to the resource grid
        reGrid(pucch3Indices) = pucch3Sym;
        reGrid(pucch3DRSIndices) = pucch3DRSSym;

        % SC-FDMA modulation
        txwave = lteSCFDMAModulate(ue, reGrid);
        
        % Channel state information: set the init time to the correct value
        % to guarantee continuity of the fading waveform after each
        % subframe
        channel.InitTime = (nsf-1)/1000;

        % Channel modeling 
        % The additional 25 samples added to the end of the waveform
        % are to cover the range of delays expected from the channel
        % modeling (a combination of implementation delay and 
        % channel delay spread)            
        rxwave = lteFadingChannel(channel, [txwave; zeros(25, ue.NTxAnts)]);           

        % Add noise at receiver
        noise = N*complex(randn(size(rxwave)), randn(size(rxwave)));
        rxwave = rxwave + noise;

        % Receiver 

        % Synchronization
        % An offset within the range of delays expected from the channel 
        % modeling (a combination of implementation delay and channel 
        % delay spread) indicates success
        offset = lteULFrameOffsetPUCCH3(ue, pucch, rxwave);
        if (offset<25)
            offsetused = offset;
        end

        % SC-FDMA demodulation
        rxgrid = lteSCFDMADemodulate(ue, rxwave(1+offsetused:end, :));

        % Channel estimation            
        [H, n0] = lteULChannelEstimatePUCCH3(ue, pucch, cec, rxgrid);
        
        % Extract REs corresponding to the PUCCH 3 from the given subframe
        % across all receive antennas and channel estimates
        [pucch3Rx, pucch3H] = lteExtractResources(pucch3Indices, rxgrid, H);

        % PUCCH 3 MMSE Equalization
        eqgrid = lteULResourceGrid(ue);    
        eqgrid(pucch3Indices) = lteEqualizeMMSE(pucch3Rx, pucch3H, n0);   

        % PUCCH 3 demodulation
        rxBits = ltePUCCH3Decode(ue, pucch, eqgrid(pucch3Indices));

        % PUCCH 3 decoding
        rxACK = lteUCI3Decode(rxBits, N_ACK);

        % Detect missed (empty rxACK) or incorrect Hybrid Automatic Repeat
        % Request (HARQ)-ACK
        % (compare against transmitted ACK)

        if (isempty(rxACK) || any(rxACK ~= ACK))
            missCount = missCount + 1;  
        end 
        
    end
    
    PMISS(nSNR) = missCount/numSubframes;

end

%% Results

plot(SNRdB, PMISS, 'b-o', 'MarkerSize', 7, 'LineWidth', 2); 
hold on;
plot(-3.7, 0.01, 'rx', 'MarkerSize', 7, 'LineWidth', 2);
xlabel('SNR (dB)');
ylabel('Probability of ACK missed detection');
title(['PUCCH format 3 ACK missed detection test' ...
        ' (TS36.104 Section 8.3.6.1)']);   
axis([SNRdB(1)-0.1 SNRdB(end)+0.1 -0.05 0.25]);
legend('simulated performance', 'target');  

%% Selected Bibliography
% # 3GPP TS36.104.

displayEndOfDemoMessage(mfilename)
