%% PDSCH 2x2 Spatial Multiplexing Throughput Conformance Test

clear;
close all;
clc;


NFrames = 200;     % Number of frames
SNRIn = 50%[2.0:1.0:25.0];    % SNR range
rng('default');  % Default random number generator seed 


% eNodeB configuration
enb = struct;                        % Create eNodeB config structure
enb.TotSubframes = 1;                % Total subframes to generate
enb.RC = 'R.11';                     % RMC number
ncw = 2;                             % Number of PDSCH codewords

% PDSCH configuration
enb.PDSCH.TxScheme = 'CDD'; %'SpatialMux';   % Transmission scheme
enb.PDSCH.RNTI = 1;                  % 16-bit UE specific mask
enb.PDSCH.Rho = -3;                  % PDSCH power adjustment factor
enb.PDSCH.CSI = 'On';                % CSI scaling of soft bits
enb.PDSCH.PMIMode = 'Wideband';      % Wideband PMI mode
enb.PDSCH.CodebookSubset = '10'; %'110000'; % Codebook subset restriction

%%
pmiDelay = 8;

%% Propagation Channel Model Configuration

cfg.Seed = 10;                   % Random channel seed
cfg.NRxAnts = 2;                 % 2 receive antennas
cfg.DelayProfile = 'EVA';        % Delay profile
cfg.DopplerFreq = 70;            % Doppler frequency
cfg.MIMOCorrelation = 'Low';     % Multi-antenna correlation
cfg.NTerms = 16;                 % Oscillators used in fading model
cfg.ModelType = 'GMEDS';         % Rayleigh fading model type
cfg.InitPhase = 'Random';        % Random initial phases
cfg.NormalizePathGains = 'On';   % Normalize delay profile power  
cfg.NormalizeTxAnts = 'On';      % Normalize for transmit antennas

%% Channel Estimator Configuration
perfectChanEstimator = false;
cec.PilotAverage = 'UserDefined';   % Type of pilot symbol averaging
cec.FreqWindow = 1;                 % Frequency window size in REs
cec.TimeWindow = 31;                % Time window size in REs
cec.InterpType = 'Cubic';           % 2D interpolation type
cec.InterpWindow = 'Centered';      % Interpolation window type
cec.InterpWinSize = 1;              % Interpolation window size

%% RMC-DL Configuration
rmc = lteRMCDL(enb, ncw);
rvSequence = rmc.PDSCH.RVSeq;
trBlkSizes = rmc.PDSCH.TrBlkSizes;

rmc.PDSCH.Modulation = {'QPSK','QPSK'};
rmc.PDSCH.CodedTrBlkSizes = [24768,26400,26400,26400,26400,0,26400,26400,26400,26400;24768,26400,26400,26400,26400,0,26400,26400,26400,26400].*0.5;

%% Setup HARQ Processes 
harqTable = hHARQTable();

%% Set Propagation Channel Model Sampling Rate
ofdmInfo = lteOFDMInfo(rmc);
cfg.SamplingRate = ofdmInfo.SamplingRate;    

%% Processing Loop
Dims = lteDLResourceGridSize(rmc);
P = Dims(3);

% Initialize variables used in the simulation and analysis
resultIndex = 1;   % Initialize SNR test point counter index
offsets = 0;       % Initialize frame offset value

% Store results for each SNR point and each subframe containing data for
% the whole simulation
nDataTBS = sum(trBlkSizes(:)~=0)*NFrames;
totalBLKCRC = zeros(numel(SNRIn), nDataTBS);
bitThroughput = zeros(numel(SNRIn), nDataTBS);

for SNRdB_index = 1:numel(SNRIn)
    fprintf('\nSimulating at %gdB SNR for %d Frame(s)\n' ,SNRIn(SNRdB_index), NFrames);

    % Initialize result store for SNR point tested
    blkCRC = zeros(nDataTBS/ncw, ncw);  % Intermediate block CRC 
    bitTput = zeros(nDataTBS/ncw, ncw); % Intermediate throughput
    dataSubframeIndex = 1;
    
    % Initialize state of all HARQ processes
%     harqprocess = hDownlinkSpatialMuxNewHARQProcess(rmc);
        harqprocess = hDownlinkNewHARQProcess(rmc);

    harqProcesses = repmat(harqprocess, 1, max(harqTable));
    
    % Use random PMIs for the first 'pmiDelay' subframes until feedback is
    % available from the UE
    pmidims = ltePMIInfo(rmc,rmc.PDSCH);
    txPMIs = randi([0 pmidims.MaxPMI], pmidims.NSubbands, pmiDelay);  

    for subframeNo = 0:(NFrames*10-1)
        
        % Update subframe number
        rmc.NSubframe = subframeNo;

        % Get HARQ index for given subframe from HARQ index table
        harqIdx = harqTable(mod(subframeNo, length(harqTable))+1);

        % Update current HARQ process
        harqProcesses(harqIdx) = hHARQScheduling( ...
            harqProcesses(harqIdx), subframeNo, rvSequence);

        % Update the PDSCH transmission config with HARQ process state
        rmc.PDSCH = harqProcesses(harqIdx).txConfig;      
        data = harqProcesses(harqIdx).data;

        % Set the PMI to the appropriate value in the delay queue
        pmiIdx = mod(subframeNo, pmiDelay);     % PMI index in delay queue
        rmc.PDSCH.PMISet = txPMIs(:, pmiIdx+1); % Set PMI

        % Create transmit waveform and add 25 sample padding
        [WAVEFORM,GRID,RMCCFGOUT] = lteRMCDLTool(rmc, data);
        figure();
        scatter(real(reshape(GRID(:,:,1),[],1)),imag(reshape(GRID(:,:,1),[],1)));
        grid on;
        
        figure();
        scatter(real(reshape(GRID(:,:,2),[],1)),imag(reshape(GRID(:,:,2),[],1)));
        grid on;
        
        txWaveform = [WAVEFORM; zeros(25, P)];

        % The initialization time for channel modeling is set each subframe
        % to simulate a continuously varying channel
        cfg.InitTime = subframeNo/1000;

        % Pass data through channel model
        rxWaveform = lteFadingChannel(cfg, txWaveform);

        % Calculate noise gain including compensation for downlink power
        % allocation
        SNR = 10^((SNRIn(SNRdB_index)-rmc.PDSCH.Rho)/20);

        % Normalize noise power to take account of sampling rate, which is
        % a function of the IFFT size used in OFDM modulation, and the 
        % number of antennas
        N0 = 1/(sqrt(2.0*rmc.CellRefP*double(ofdmInfo.Nfft))*SNR);

        % Create additive white Gaussian noise
        noise = N0*complex(randn(size(rxWaveform)), ...
                            randn(size(rxWaveform)));

        % Add AWGN to the received time domain waveform        
        rxWaveform = rxWaveform + noise;

        % Once every frame, on subframe 0, calculate a new synchronization
        % offset
        if (mod(subframeNo,10) == 0)
            offset = lteDLFrameOffset(rmc, rxWaveform);
            if (offset > 25)
                offset = offsets(end);
            end
            offsets = [offsets offset]; %#ok
        end
        
        % Synchronize the received waveform
        rxWaveform = rxWaveform(1+offset:end, :);

        % Perform OFDM demodulation on the received data to recreate the
        % resource grid
        rxSubframe = lteOFDMDemodulate(rmc, rxWaveform);

        % Channel estimation
        if(perfectChanEstimator)
            estChannelGrid = lteDLPerfectChannelEstimate(rmc, cfg, offset); %#ok
            n = lteOFDMDemodulate(rmc, noise(1+offset:end ,:));
            noiseEst = var(n(:));
        else
            [estChannelGrid, noiseEst] = lteDLChannelEstimate( ...
                rmc, cec, rxSubframe);
        end

        % Decode the PDSCH. The received subframe is scaled to counteract
        % any attenuation applied at the transmitter.
        rxEncodedBits = ltePDSCHDecode(rmc, rmc.PDSCH,...
              rxSubframe*(10^(-rmc.PDSCH.Rho/20)), ...
              estChannelGrid, noiseEst);

        % Extract the current subframe transport block size(s)
        TBSs = trBlkSizes(:, mod(subframeNo, 10)+1).';

        % Decode the DL-SCH
        [decbits, harqProcesses(harqIdx).blkerr, ...
            harqProcesses(harqIdx).decState] = lteDLSCHDecode( ...
            rmc, rmc.PDSCH, TBSs, rxEncodedBits, ...
            harqProcesses(harqIdx).decState);
        
        % Store block CRC and throughput results for subframes containing
        % transport data, ignore subframes containing no data
        if(any(TBSs ~= 0))
            blkCRC(dataSubframeIndex,:) = harqProcesses(harqIdx).blkerr;
            bitTput(dataSubframeIndex,:) = ...
                TBSs.*(1-harqProcesses(harqIdx).blkerr);
            dataSubframeIndex = dataSubframeIndex + 1;
        end

        % Provide PMI feedback to the eNodeB
        PMI = ltePMISelect(rmc, rmc.PDSCH, estChannelGrid, noiseEst);
        txPMIs(:, pmiIdx+1) = PMI;

    end

    % Record the block CRC error and bit throughput for the total number of
    % frames simulated at an SNR point
    totalBLKCRC(resultIndex,:) = blkCRC(:);
    bitThroughput(resultIndex,:) = bitTput(:);
    resultIndex = resultIndex + 1;

    % Display the results dynamically in the command window
    fprintf([['Throughput(Mbps) for ', num2str(NFrames) ' Frame(s) '],...
              '= %.4f\n'], mean(bitTput(:))* (nDataTBS/(NFrames*10))/1e3);
    fprintf(['Throughput(%%) for ', num2str(NFrames) ' Frame(s) = %.4f\n'],...
             100*(1-mean(blkCRC(:))));
end

%% Throughput Results
% The throughput results for the simulation are displayed in the MATLAB(R)
% command window after each SNR point is completed, and are also captured
% in |totalBLKCRC| and |bitThroughput|. |totalBLKCRC| is a matrix where
% each row contains the results of decoding the CRC for a particular SNR.
% Each column contains the CRC result for a transport block containing
% PDSCH data at an SNR. |bitThroughput| is a matrix where each row contains
% the bit throughput per subframe for a particular SNR. Each column
% contains the throughput result for a transport block containing PDSCH
% data at an SNR.

%% Results


residualBler = mean(totalBLKCRC,2);
figure();
semilogy(SNRIn,residualBler,'r-o');
 title([ num2str(length(rvSequence))  'Tx residualBler for ', num2str(NFrames) ' Frames']);
    xlabel('SNR (dB)'); ylabel('residualBler');
grid on;
set(gca,  'XLim', [min(SNRIn) max(SNRIn)],'YLim', [0.0001 1]);



figure;
plot(SNRIn,mean(bitThroughput,2)*0.9*2,'-*');   % adjust ,(12960*3+9528*2)*3 / (12960*7+9528*8) = 1.0411, *0.5*1.0411
% axis([-5 3 200 400])
title(['Throughput for ', num2str(NFrames) ' frame(s)'] );
xlabel('SNRdB'); ylabel('Throughput (kbps)');
grid on;
hold on;
plot(SNRIn,mean(mean(trBlkSizes,2),1)*0.7*2*ones ...
    (1,numel(SNRIn)),'--rs');
legend('Simulation Result','70 Percent Throughput','Location','SouthEast');

% Second graph shows the total throughput as a percentage of CRC passes against SNR range
figure;
plot(SNRIn,100*(1-mean(totalBLKCRC,2)),'-*');
% axis([-5 3 50 110])
title(['Throughput for ', num2str(NFrames) ' frame(s)'] );
xlabel('SNRdB'); ylabel('Throughput (%)');
grid on;
hold on;
plot(SNRIn,70*ones(1,numel(SNRIn)),'--rs');
legend('Simulation Result','70 Percent Throughput','Location','SouthEast');

close all;
