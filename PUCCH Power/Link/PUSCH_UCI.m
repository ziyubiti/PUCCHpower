%% PUSCH HARQ-ACK Detection Conformance Test
% * _ACK false detection 
% * _The ACK missed detection 
% ack is transmitted on odd subframes and doing missed dectecion
% in even subframe ,doing false detection;
% 

%% Simulation Configuration
clc;
close all;
clear;
path(path,'../Model');
numSubframes = 2000;  % Number of frames to simulate at each SNR
SNRdB = [-1.0:1.0: 10.0  ]; % SNR points to simulate

%% UE Configuration
frc.TotSubframes = 1;   % Total number of subframes to generate
frc.NCellID = 10;       % Cell identity
frc.RC = 'A3-1';        % FRC number
frc = lteRMCUL(frc);
frc.NULRB = 100;
frc.PUSCH.BetaACK = 12.625;  % etu setting 5(6.25), eva setting 8(12.625)

%% Propagation Channel Model Configuration

chcfg.NRxAnts = 2;               % Number of receive antennas
chcfg.DelayProfile = 'ETU';      % Delay profile
chcfg.DopplerFreq = 70;          % Doppler frequency    
chcfg.MIMOCorrelation = 'Low';   % MIMO correlation
chcfg.Seed = 65535;                 % Channel seed
chcfg.NTerms = 16;               % Oscillators used in fading model
chcfg.ModelType = 'GMEDS';       % Rayleigh fading model type
chcfg.InitPhase = 'Random';      % Random initial phases
chcfg.NormalizePathGains = 'On'; % Normalize delay profile power 
chcfg.NormalizeTxAnts = 'On';    % Normalize for transmit antennas

% Set channel model sampling rate
info = lteSCFDMAInfo(frc);
chcfg.SamplingRate = info.SamplingRate;

%% Channel Estimator Configuration
cec.PilotAverage = 'UserDefined'; % Type of pilot averaging
cec.FreqWindow = 9;               % Frequency averaging windows in REs
cec.TimeWindow = 1;               % Time averaging windows in REs
cec.InterpType = 'cubic';         % Interpolation type
cec.Reference = 'Antennas';       % Reference for channel estimation 

%% Loop for SNR Values
% Initialize variables used in the simulation and analysis
pFalse = zeros(size(SNRdB));  % Probability of false detection at each SNR
pMissed = zeros(size(SNRdB)); % Probability of missed detection at each SNR

for nSNR = 1:length(SNRdB)

    % Initialize the random number generator stream
    rng('default'); 
    
    fprintf('\nSimulating at %g dB SNR for a total %d SubFrame(s)', ...
        SNRdB(nSNR), numSubframes);   
    
    % Extract SNR to test
    SNR = 10^(SNRdB(nSNR)/20);
    
    % Scale noise to ensure the desired SNR after SC-FDMA demodulation
    N = 1/(SNR*sqrt(double(info.Nfft)))/sqrt(2.0);    

    offsetUsed = 0;
    falseCount = 0;
    falseTests = 0;
    missCount = 0;
    missTests = 0;

    for subframeNo = 0:(numSubframes-1)

        % Updating subframe number
        frc.NSubframe = mod(subframeNo, 10);

        % Transmit ACK on every odd subframe
        if (mod(subframeNo, 2)==0)
            ACK = [];
            falseTests = falseTests + 1;
        else
            ACK = 1;
            missTests = missTests + 1;
        end            

        % Create random data to transmit
        trblklen = frc.PUSCH.TrBlkSizes(frc.NSubframe+1);
        trblk = randi([0 1], trblklen, 1);                                               

        % [WAVEFORM,GRID,RMCCFGOUT] = lteRMCULTool(RMCCFG,TRDATA,CQI,RI,ACK)
        txWaveform = [lteRMCULTool(frc, trblk, [], [], ACK); zeros(25, 1)];        

        % Pass waveform through fading channel model
        chcfg.InitTime = subframeNo/1000;
        rxWaveform = lteFadingChannel(chcfg, txWaveform);

        % Add noise
        noise = N*complex(randn(size(rxWaveform)), ...
            randn(size(rxWaveform)));
        rxWaveform = rxWaveform + noise;

        % Synchronization
        offset = lteULFrameOffset(frc, frc.PUSCH,rxWaveform);
        if (offset < 25)
            offsetUsed = offset;
        end

        % SC-FDMA demodulation
        rxSubframe = lteSCFDMADemodulate(frc, ...
            rxWaveform(1+offsetUsed:end, :));

        % Channel Estimation
        [estChannelGrid, noiseEst] = lteULChannelEstimate( ...
            frc, frc.PUSCH, cec, rxSubframe);

        % PUSCH indices for given subframe
        puschIndices = ltePUSCHIndices(frc,frc.PUSCH);

        % Minimum Mean Squared Error (MMSE) equalization
        rxSymbols = lteEqualizeMMSE(rxSubframe(puschIndices), ...
                        estChannelGrid(puschIndices), noiseEst);
     
        frc.PUSCH = lteULSCHInfo(frc, frc.PUSCH, trblklen, ...
            0, 0, 1, 'chsconcat');
  
        rxEncodedBits = ltePUSCHDecode(frc, frc.PUSCH, rxSymbols); 

        % UL-SCH channel deinterleaving
        [deinterleavedData, ccqi, cri, cack] = ...
            lteULSCHDeinterleave(frc, frc.PUSCH, rxEncodedBits);

        % HARQ-ACK decoding
        rxACK = lteACKDecode(frc.PUSCH, cack);
        
        % Detect false or missed HARQ-ACK
        if (isempty(ACK) && ~isempty(rxACK))
            falseCount = falseCount + 1;
        end            
        if (~isempty(ACK) && ~isequal(ACK,rxACK)) 
            missCount = missCount + 1;
        end       
          
    end     

    % Calculate false or missed HARQ-ACK probability
    pFalse(nSNR) = falseCount/falseTests;
    pMissed(nSNR) = missCount/missTests;

end

%% Analysis
hHARQACKResultsV2(SNRdB, pFalse, pMissed);

