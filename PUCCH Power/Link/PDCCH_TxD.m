%% PDCCH Conformance Test

clear;
clc;
close all;

nFrames = 10000;                      % Number of frames to simulate
snrIn = [-6.0:1.0:-1.0]; % SNR points to simulate

%% PDCCH/PCFICH Configuration

rmc = lteRMCDL('R.11','TDD');

rmc.TotSubframes = 1;

% Setup DCI and PDCCH according to TS 36.101 Section 8.4.1.1 for R.15 FDD
% RMC which requires the CFI, DCI format, HICH group multiplier,
% aggregation level and power to be set as defined for the test
rmc.CFI = 2;                     % OFDM symbols for PDCCH
rmc.Ng = 'One';                  % HICH group multiplier
rmc.PDSCH.DCIFormat = 'Format2'; % Set the DCI Format
rmc.PDSCH.PDCCHFormat = 2;       % Set the aggregation level to be 8
rmc.PDSCH.PDCCHPower = -3;        % Relative power is 0dB for single port 
rmc.CellRefP = 2;                   % Number of transmit antenna ports

rmc.PDSCH.TxScheme = 'TxDiversity'; % PDSCH transmission scheme
rmc.PDSCH.NLayers = 2;              % Number of transmission layers

% Setup the OCNG to fill all unused data and control region REs
rmc.OCNGPDSCH.TxScheme = 'TxDiversity'; % PDSCH OCNG transmission scheme
rmc.OCNGPDSCHEnable = 'On';
rmc.OCNGPDSCHPower = -3;
rmc.OCNGPDCCHEnable = 'On';
rmc.OCNGPDCCHPower = -3;

%% Propagation Channel Model Configuration
cfg = struct;                    % Initialize channel config structure
cfg.Seed = 4;                    % Random channel seed
cfg.NRxAnts = 2;                 % 2 receive antennas
cfg.DelayProfile ='EVA';         % Delay profile
cfg.DopplerFreq = 70;            % Doppler frequency in Hz
cfg.MIMOCorrelation = 'Low';     % Multi-antenna correlation
cfg.NTerms = 16;                 % Oscillators used in fading model
cfg.ModelType = 'GMEDS';         % Rayleigh fading model type
cfg.InitPhase = 'Random';        % Random initial phases
cfg.NormalizePathGains = 'On';   % Normalize delay profile power     
cfg.NormalizeTxAnts = 'On';      % Normalize for transmit antennas

%% Channel Estimator Configuration
perfectChannelEstimator = false;

% Configure channel estimator
cec.PilotAverage = 'UserDefined';   % Type of pilot symbol averaging
cec.FreqWindow = 1;                 % Frequency window size in REs
cec.TimeWindow = 31;                % Time window size in REs
cec.InterpType = 'Cubic';           % 2D interpolation type
cec.InterpWindow = 'Centered';      % Interpolation window type
cec.InterpWinSize = 3;              % Interpolation window size


%% Processing 

% dims is a three element vector [K L P]: where K = no of subcarriers,
% L = no of OFDM symbols and P = no of transmit antennas
dims = lteDLResourceGridSize(rmc);  % Determine resource grid dimensions
L = dims(2);                        % Number of OFDM symbols per subframe
P = dims(3);                        % Number of transmit antennas

% Initialize the variables used in the simulation and analysis
resultIndex = 1;                    % Initialize SNR index
initOffsets = 0;                    % Initialize frame offset 

% Set the random number generator to default value
rng('default');

% Number of subframes carrying transport blocks with data and control
% within a frame. For RMC R.15 this should be 9, as out of the 10 subframes
% in a frame subframe 5 does not carry any data or control
% nDataTBSPerFrame = sum(rmc.PDSCH.TrBlkSizes(:) ~= 0);
  nDataTBSPerFrame = 1;  % only in sf0

% Total Pmdsg vector
totalPmdsg = zeros(numel(snrIn),nDataTBSPerFrame*nFrames); 

for snrdb = snrIn
    fprintf('\nSimulating at %gdB SNR for a total %d Frame(s)\n',...
    snrdb,nFrames);

    % Initialize result vector for each SNR point
    totalPmdsgSNR = zeros(nFrames,nDataTBSPerFrame);% Total block CRC 
    
    % Initialize offset vector for each frame
    offsets = initOffsets;
    
    for FrameNo = 1:nFrames
        
        % Set the subframe number
        rmc.NSubframe = 0;
        
        % Generate the waveform
        [txWaveform,txGrid,info] = lteRMCDLTool(rmc,[1 0 0 1]);
        
        % Initialize result store for frame
        pmdsg = zeros(nDataTBSPerFrame,1);   
        dataSubframeIndex = 1;

        % Set sampling rate of channel to that of OFDM modulation     
        cfg.SamplingRate = info.SamplingRate;

        % Set channel offset to current frame (1 frame = 10ms)
        cfg.InitTime = (FrameNo-1)*(rmc.TotSubframes)/1000;

        % Pass data through the fading channel model. 
        % An additional 25 samples are added to the end of the waveform.
        % These are to cover the range of delays expected from the channel
        % modeling (a combination of implementation delay and channel
        % delay spread).
        rxWaveform = lteFadingChannel(cfg,[txWaveform ; zeros(25,P)]);

        % Noise setup
        SNR = 10^((snrdb-rmc.PDSCH.PDCCHPower)/20);    % Linear SNR

        % Normalize noise power to take account of sampling rate, which is
        % a function of the IFFT size used in OFDM modulation, and the 
        % number of antennas
        N0 = 1/(sqrt(2.0*rmc.CellRefP*double(info.Nfft))*SNR); 

        % Create additive white Gaussian noise
        noise = N0*complex(randn(size(rxWaveform)), ...
                            randn(size(rxWaveform)));

        % Add AWGN to the received time domain waveform        
        rxWaveform = rxWaveform + noise;

        % Perform receiver synchronization 
        offset = lteDLFrameOffset(rmc,rxWaveform);

        % Determine if frame synchronization was successful
        if (offset > 25)
            offset = offsets(end);
        else
            offsets = [offsets offset]; %#ok
        end
        if (offset>0)
            rxWaveform = rxWaveform(1+offset:end,:);
        end

        % Perform OFDM demodulation on the received data to recreate the
        % resource grid
        rxGrid = lteOFDMDemodulate(rmc,rxWaveform);
        
        % Channel estimation
        if(perfectChannelEstimator)
            estChannelGrid = lteDLPerfectChannelEstimate(rmc, cfg, offset);  %#ok<UNRCH>
            n = lteOFDMDemodulate(rmc, noise(1+offset:end ,:));
            noiseest = var(n(:));
        else
            [estChannelGrid, noiseest] = lteDLChannelEstimate( ...
                rmc, cec, rxGrid);
        end

        % Process subframes 0 to 9 within received frame. Test for PCFICH
        % and if successful, test for PDCCH. If both PCFICH and PDCCH are
        % successfully received, the transmission is successful
        for sf = 0:rmc.TotSubframes-1;

            % Increment subframe number for correct decoding of data
            rmc.NSubframe = mod(sf,10);
            
            % No data is transmitted in subframe 5, so skip this
            if(rmc.NSubframe ~= 5) 
            
                % Extract one subframe for analyzing at a time from the
                % received grid        
                rxSubframe = rxGrid(:,L*sf+1:L*(sf+1),:);

                % Extract the estimated channel response for the subframe
                % being analyzed
                chSubframe = estChannelGrid(:,L*sf+1:L*(sf+1),:,:);

                % Perform PCFICH decoding of received data using the estimate
                % of the channel
                pcfichIndices = ltePCFICHIndices(rmc);
                [rxPcfichSym,pcfichHestSym] = lteExtractResources(pcfichIndices,rxSubframe,chSubframe);
                rxPcfich = ltePCFICHDecode(rmc,rxPcfichSym,pcfichHestSym,noiseest);
                rxCFI = lteCFIDecode(rxPcfich);

                if rxCFI ~= rmc.CFI
                    rxFail = 1;
                else
                   % CFI decoded fine, now check if DCI can be decoded

                   % Extract and decode PDCCH bits
                   pdcchIndices = ltePDCCHIndices(rmc);
                   [rxPdcchSym,pdcchHestSym] = lteExtractResources(pdcchIndices,rxSubframe,chSubframe);
                   rxPdcchBits = ltePDCCHDecode(rmc,rxPdcchSym,pdcchHestSym,noiseest); 

                   % PDCCH blind search, demask PDCCH candidate using RNTI
                    ueConfig.RNTI = rmc.PDSCH.RNTI; 
                   [rxDCI,rxDCIBits] = ltePDCCHSearch(rmc,ueConfig,rxPdcchBits);

                   if ~isempty(rxDCI) && isfield(rxDCI{1},'DCIFormat') && isequal(rxDCI{1}.DCIFormat,rmc.PDSCH.DCIFormat)
                       rxFail = 0;
                   else
                       rxFail = 1;
                   end
            
                end
                % Store the values and increment the subframe index
                pmdsg(dataSubframeIndex, :) = rxFail;
                dataSubframeIndex = dataSubframeIndex + 1;
            end
        end
        % Store resulting pmdsg results 
        totalPmdsgSNR(FrameNo,:) = pmdsg(:);

    end
    % Record the missed detection for the total number of frames simulated
    % at a particular SNR
    totalPmdsg(resultIndex,:) = totalPmdsgSNR(:); % CRC 
    resultIndex = resultIndex + 1;
    
end

%% Plot Results

% Pm-dsg should be less than 1% at -1.7 dB SNR for single port
pmdsgVal = mean(totalPmdsg,2);
figure();
semilogy(snrIn,pmdsgVal,'-*');
title(['Pm-dsg for ', num2str(nFrames) ' frames'] );
xlabel('SNR(dB)'); ylabel('Pm-dsg');
grid on;
legend('PDCCH format2','Location','NorthEast');
% ylim([fix(min(pmdsgVal)-2) fix(max(pmdsgVal)+2)]);

