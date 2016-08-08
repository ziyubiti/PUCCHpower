%%  PDSCH Transmit Diversity Throughput Conformance Test
clear;
clc;
close all;
path(path,'../Model');

NFrames = 1000;                        % Number of frames
SNRdB = [-10:4:6];         % SNR range

% eNodeB Configuration
enb = struct;                        % eNodeB config structure
enb.TotSubframes = 1;                % Total subframes RMC will generate
enb.RC = 'R.11';                     % RMC number
enb.DuplexMode = 'FDD';
enb.TDDConfig = 1;
enb.SSC = 4;

% Channel Configuration
channel = struct;                    % Channel config structure
channel.Seed = 2;                    % Random channel seed
channel.NRxAnts = 2;                 % 2 receive antennas
channel.DelayProfile ='EPA';         % Delay profile
channel.DopplerFreq = 5;             % Doppler frequency
channel.MIMOCorrelation = 'Medium';  % Multi-antenna correlation
channel.NTerms = 16;                 % Oscillators used in fading model
channel.ModelType = 'GMEDS';         % Rayleigh fading model type
channel.InitPhase = 'Random';        % Random initial phases
channel.NormalizePathGains = 'On';   % Normalize delay profile power 
channel.NormalizeTxAnts = 'On';      % Normalize for transmit antennas

% Channel Estimator Configuration
cec = struct;                        % Channel estimation config structure
cec.PilotAverage = 'UserDefined';    % Type of pilot symbol averaging
cec.FreqWindow = 15;                 % Frequency window size
cec.TimeWindow = 15;                 % Time window size
cec.InterpType = 'Cubic';            % 2D interpolation type
cec.InterpWindow = 'Centered';       % Interpolation window type
cec.InterpWinSize = 1;               % Interpolation window size

% PDSCH Configuration
enb.PDSCH.TxScheme = 'TxDiversity';  % Transmission scheme
enb.PDSCH.RNTI = 1;                  % 16-bit User Equipment (UE) mask
enb.PDSCH.Rho = -3;                  % PDSCH RE power adjustment factor
enb.PDSCH.CSI = 'Off';               % No CSI scaling of soft bits

% Simulation Variables
totalBLKCRC = [];                    % Define total block CRC vector
bitThroughput = [];                  % Define total bit throughput vector

%% System Processing


% Generate the RMC configuration structure for RMC R.12
rmc = lteRMCDL(enb);
rvSeq = 0;%[ 0 2 3 1];%[0 2 3 1];
rmc.PDSCH.RVSeq = rvSeq;

rmc.NDLRB = 100;
rmc.PDSCH.Modulation = 'QPSK';
rmc.PDSCH.PRBSet = [0:1;0:1]';
rmc.PDSCH.TrBlkSizes = [208,208,208,208,208,0,208,208,208,208];
rmc.PDSCH.CodedTrBlkSizes = [528,528,528,528,528,0,528,528,528,528];

    
% Transport block sizes for each subframe in a frame
trBlkSizes = rmc.PDSCH.TrBlkSizes;
codedTrBlkSizes = rmc.PDSCH.CodedTrBlkSizes;
    
% Determine resource grid dimensions
dims = lteDLResourceGridSize(rmc);
p = dims(3);
    
% Set up channel model sampling rate
ofdmInfo = lteOFDMInfo(rmc);
channel.SamplingRate = ofdmInfo.SamplingRate;    

% Generation HARQ table for 8-HARQ processes
harqTable = hHARQTableV2(enb.DuplexMode);

% Initializing state of all HARQ processes
for i=1:9    % 9 for FDD
     harqProcess_init(i) = hTxDiversityNewHARQProcess ...
                    (trBlkSizes(i),codedTrBlkSizes(i),rvSeq); %#ok<SAGROW>
end

% The temporary variables 'rmc_init' and 'channel_init' are used to create
% the temporary variables 'rmc' and 'channel' within the SNR loop to create
% independent simulation loops for the parfor loop
rmc_init = rmc;
channel_init = channel;

% bj ; add multi-user and multi-cell
usersPDSCHpower = [0 6 -3 3];
usersNCellID = [3 150 150 3] ;
usersRNTI = [95 3  202 37] ;
% userCyclicShift = [0 3 6 1];
% userSeqGroup = [0 1  2 3];
ueChannelSeed = [ 2 7 200 500];




% Main loop
fprintf('\nSimulating %d SNR points for %d frame(s) each\n', ...
            length(SNRdB),NFrames);
        
for index = 1:numel(SNRdB)
% To enable the use of parallel computing for increased speed comment out
% the 'for' statement above and uncomment the 'parfor' statement below.
% This needs the Parallel Computing Toolbox. If this is not installed
% 'parfor' will default to the normal 'for' instruction.
%parfor index = 1:numel(SNRdB)

    
    % Set the random number generator seed depending to the loop variable
    % to ensure independent random streams
    rng(index,'combRecursive');
    fprintf('\nSimulating at %g dB SNR for total %d Frames', SNRdB(index),NFrames);

    % Set up variables for the SNR loop
    offsets = 0;    % Initialize overall frame offset value for the SNR
    offset = 0;     % Initialize frame offset value for the radio frame
    rmc = rmc_init; % Initialize RMC configuration
    channel = channel_init; % Initialize channel configuration 
    blkCRC = [];    % Define intermediate block CRC vector
    bitTput = [];   % Intermediate bit throughput vector
              
    % Initializing state of all HARQ processes
    harqProcess = harqProcess_init;
    
    for subframeNo = 0:(NFrames*10-1)

        % Updating subframe number
        rmc.NSubframe = subframeNo;
        
        for user = 1:1
            rmc.NCellID = usersNCellID(user);            
            rmc.PDSCH.RNTI = usersRNTI(user);
            

        % HARQ index table
        harqIdx = harqTable(mod(subframeNo,length(harqTable))+1); 

        % Update HARQ process
        harqProcess(harqIdx) = hTxDiversityHARQScheduling( ...
                               harqProcess(harqIdx)); 

        % Updating the RV value for correct waveform generation
        rmc.PDSCH.RV = harqProcess(harqIdx).rvSeq ...
                       (harqProcess(harqIdx).rvIdx);
                   
        rmc.PDSCH.RVSeq = harqProcess(harqIdx).rvSeq ...
                          (harqProcess(harqIdx).rvIdx);
                      
        txWaveform = lteRMCDLTool(rmc, ...
                   {harqProcess(harqIdx).dlschTransportBlk});
%         [txWaveform, txSubframe] = lteRMCULTool( ...
%                 frc, harqProc(harqIdx).ulschTransportBlk);
               
        % Initialize at time zero  
        channel.InitTime = subframeNo/1000;

        % Pass data through the fading channel model 
        if (user==1)
            rxWaveform = lteFadingChannel(channel,[txWaveform;zeros(25,p)]);
        else
            rxWaveform = rxWaveform + ...
                lteFadingChannel(channel,[txWaveform*10^(usersPDSCHpower(user)/20); zeros(25,p)]);
        end;
        
        end;   % end user

        % Noise setup including compensation for downlink power allocation
        SNR = 10^((SNRdB(index)-rmc.PDSCH.Rho)/20);    % Linear SNR

        % Normalize noise power to take account of sampling rate, which is
        % a function of the IFFT size used in OFDM modulation, and the 
        % number of antennas
        N0 = 1/(sqrt(2.0*rmc.CellRefP*double(ofdmInfo.Nfft))*SNR); 

        % Create additive white Gaussian noise
        noise = N0*complex(randn(size(rxWaveform)), ...
                            randn(size(rxWaveform)));
                        
        % Add AWGN to the received time domain waveform
        rxWaveform = rxWaveform + noise;
        
                

        % Receiver
        % Perform synchronization
        % An offset within the range of delays expected from the channel 
        % modeling(a combination of implementation delay and channel delay
        % spread) indicates success
        if (mod(subframeNo,10)==0)
            [offset] = lteDLFrameOffset(rmc,rxWaveform);
            if (offset > 25)
                offset = offsets(end);
            end
            offsets = [offsets offset];  %#ok<AGROW>
        end
        rxWaveform = rxWaveform(1+offset:end,:);                        

        % Perform OFDM demodulation on the received data to recreate the
        % resource grid
        rxSubframe = lteOFDMDemodulate(rmc,rxWaveform);

        % Equalization and channel estimation
        [estChannelGrid,noiseEst] = lteDLChannelEstimate(rmc,cec, ...
                                                         rxSubframe);

        % Perform deprecoding, layer demapping, demodulation and 
        % descrambling on the received data using the estimate of
        % the channel
        rxEncodedBits = ltePDSCHDecode(rmc,rmc.PDSCH,rxSubframe* ...
                        (10^(-rmc.PDSCH.Rho/20)),estChannelGrid,noiseEst);

        % Decode DownLink Shared Channel (DL-SCH)
        [decbits,harqProcess(harqIdx).crc,harqProcess(harqIdx).decState] = ...  
            lteDLSCHDecode(rmc,rmc.PDSCH,harqProcess(harqIdx).trBlkSize, ...
            rxEncodedBits{1},harqProcess(harqIdx).decState);
    
        if(harqProcess(harqIdx).trBlkSize ~= 0)
            blkCRC = [blkCRC harqProcess(harqIdx).crc]; %#ok<AGROW>
            bitTput = [bitTput harqProcess(harqIdx).trBlkSize.*(1- ...
                      harqProcess(harqIdx).crc)]; %#ok<AGROW>
        end
    end
    % Record the block CRC and bit throughput for the total number of
    % frames simulated at a particular SNR
    totalBLKCRC(index,:) = blkCRC; %#ok<SAGROW>                    % TDD FRC, 12960 should be 3 in a Frame, but it is 2 3 2 in 3Frame,
    bitThroughput(index,:) = bitTput; %#ok<SAGROW>

end


%% Results


residualBler = mean(totalBLKCRC,2);
figure();
semilogy(SNRdB,residualBler,'r-o');
 title([ num2str(length(rvSeq))  'Tx residualBler for ', num2str(NFrames) ' Frames']);
    xlabel('SNR (dB)'); ylabel('residualBler');
grid on;
set(gca,  'XLim', [min(SNRdB) max(SNRdB)],'YLim', [0.0001 1]);



figure;
plot(SNRdB,mean(bitThroughput,2)*0.9,'-*');   % adjust ,(12960*3+9528*2)*3 / (12960*7+9528*8) = 1.0411, *0.5*1.0411
% axis([-5 3 200 400])
title(['Throughput for ', num2str(NFrames) ' frame(s)'] );
xlabel('SNRdB'); ylabel('Throughput (kbps)');
grid on;
hold on;
plot(SNRdB,mean([trBlkSizes(1:10)])*0.7*ones ...
    (1,numel(SNRdB)),'--rs');
legend('Simulation Result','70 Percent Throughput','Location','SouthEast');

% Second graph shows the total throughput as a percentage of CRC passes 
% against SNR range
figure;
plot(SNRdB,100*(1-mean(totalBLKCRC,2)),'-*');
% axis([-5 3 50 110])
title(['Throughput for ', num2str(NFrames) ' frame(s)'] );
xlabel('SNRdB'); ylabel('Throughput (%)');
grid on;
hold on;
plot(SNRdB,70*ones(1,numel(SNRdB)),'--rs');
legend('Simulation Result','70 Percent Throughput','Location','SouthEast');

close all;
