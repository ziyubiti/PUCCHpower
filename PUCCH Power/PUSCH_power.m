%% PUSCH 
clear;
close all;
clc;

numSubframes = 100;                % Number of frames to simulate at each SNR
% SNRIn = [-4.1, -2.0, 0.1];  % SNR points to simulate
% SNRIn = [-6.0, -5.0, -4.0,-3.0,-2.0,-1.0, 0.0];  % SNR points to simulate
SNRIn = 10;%[-10.0:1.0:3.0];

ue.TotSubframes = 1; % Total number of subframes to generate a waveform for
ue.NCellID = 10;     % Cell identity
% ue.RC = 'A3-2';      % FRC number
ue.RC = 'A3-7';      % FRC number

%% Propagation Channel Model Configuration

chcfg.NRxAnts = 1;               % Number of receive antenna
chcfg.DelayProfile = 'EPA';      % Delay profile
chcfg.DopplerFreq = 5.0;         % Doppler frequency    
chcfg.MIMOCorrelation = 'Low';   % MIMO correlation
chcfg.Seed = 100;                % Channel seed    
chcfg.NTerms = 16;               % Oscillators used in fading model
chcfg.ModelType = 'GMEDS';       % Rayleigh fading model type 
chcfg.InitPhase = 'Random';      % Random initial phases
chcfg.NormalizePathGains = 'On'; % Normalize delay profile power
chcfg.NormalizeTxAnts = 'On';    % Normalize for transmit antennas

%% Channel Estimator Configuration
% Channel estimation settings are defined using a structure.

cec.FreqWindow = 13;              % Frequency averaging windows in REs
cec.TimeWindow = 1;               % Time averaging windows in REs
cec.InterpType = 'cubic';         % Interpolation type
cec.PilotAverage = 'UserDefined'; % Type of pilot averaging 
cec.Reference = 'Antennas';       % Reference for channel estimation

%% Uplink RMC Configuration
% Generate FRC configuration structure for A3-2
frc = lteRMCUL(ue);

rvSeq = [0];
frc.PUSCH.RVSeq = rvSeq;

% bj: cofig para
% frc.DuplexMode = 'TDD';
 frc.PUSCH.PRBSet = [10:19;10:19].';    % 0-99 prb index
 frc.PUSCH.TrBlkSizes = [872,872,872,872,872,872,872,120,872,872];  %10RB
%   frc.PUSCH.TrBlkSizes = [328,328,328,328,328,328,328,328,328,328];   % 4RB
%    frc.PUSCH.TrBlkSizes = [72,72,72,72,72,72,72,72,72,72];   % 1RB
 frc.PUSCH.CodedTrBlkSizes = [2880,2880,2880,2880,2880,2880,2880,2880,2880,2880]; % 10RB
% frc.PUSCH.CodedTrBlkSizes = [1152,1152,1152,1152,1152,1152,1152,1152,1152,1152];  % 4RB
% frc.PUSCH.CodedTrBlkSizes = [288,288,288,288,288,288,288,288,288,288];  % 1RB

% Transport block sizes for each subframe within a frame
trBlkSizes = frc.PUSCH.TrBlkSizes;
codedTrBlkSizes = frc.PUSCH.CodedTrBlkSizes;




%% Setup HARQ Processes 
% Generate HARQ process table
noHarqProcesses = 8;
harqTable = mod(0:noHarqProcesses-1, noHarqProcesses)+1;  

%% Set Propagation Channel Model Sampling Rate
info = lteSCFDMAInfo(frc);
chcfg.SamplingRate = info.SamplingRate;     


%% Processing Loop
% Initialize variables used in the simulation and analysis
totalBLKCRC = zeros(numel(SNRIn), numSubframes);   % Total block CRC vector
bitThroughput = zeros(numel(SNRIn), numSubframes); % Total throughput vector
resultIndex = 1;        % Initialize frame counter index

% bj 
usersPUSCHpower = [0 6 -3 3];
usersNCellID = [3 150 150 3] ;
usersRNTI = [95 3  202 37] ;
userCyclicShift = [0 3 6 1];
userSeqGroup = [0 1  2 3];
ueChannelSeed = [ 2 7 200 500];

Sp = zeros(numSubframes,1);
Sp0 = zeros(numSubframes,1);
Spave = zeros(size(SNRIn));
Spstd = zeros(size(SNRIn));
Sp0ave = zeros(size(SNRIn));
Sp0std = zeros(size(SNRIn));

for nSNR = 1:length(SNRIn) %SNRdB = SNRIn

   
    
    % Calculate required AWGN channel noise
%     SNR = 10^(SNRdB/20);
    SNR = 10^(SNRIn(nSNR)/20);     
    N = 1/(SNR*sqrt(double(info.Nfft)))/sqrt(2.0);    
    rng('default');
    
     fprintf('\nSimulating at %g dB SNR for a total %d SubFrame(s)', ...
        SNRIn(nSNR), numSubframes);
    
    % Store results for every subframe at SNR point
    bitTp = zeros(1,numSubframes);  % Intermediate bit throughput vector	
    blkCRC = zeros(1, numSubframes); % Intermediate block CRC vector         
    
    % Initialize state of all HARQ processes
    for i = 1:8
        harqProc(i) = hPUSCHNewHARQProcess( ...
            trBlkSizes(i), codedTrBlkSizes(i), rvSeq); %#ok
    end

    offsetused = 0;
    for subframeNo = 7%0:numSubframes-1%(NFrames*10-1)

        % Update subframe number
        frc.NSubframe = subframeNo;
        
        for user = 1:1
            frc.NCellID = usersNCellID(user);
            frc.CyclicShift = userCyclicShift(user);
            frc.RNTI = usersRNTI(user);
            frc.SeqGroup = userSeqGroup(user);
            
            
            % Get HARQ index for given subframe from HARQ index table
            harqIdx = harqTable(mod(subframeNo, length(harqTable))+1);
            
            % Update current HARQ process
            harqProc(harqIdx) = hPUSCHHARQScheduling(harqProc(harqIdx));
            frc.PUSCH.RV = harqProc(harqIdx).rvSeq(harqProc(harqIdx).rvIdx);
            frc.PUSCH.RVSeq = harqProc(harqIdx).rvSeq(harqProc(harqIdx).rvIdx);
            
            % Create an SC-FDMA modulated waveform
            [txWaveform, txSubframe] = lteRMCULTool( ...
                frc, harqProc(harqIdx).ulschTransportBlk);
            
            % Transmit an additional 25 samples at the end of the waveform to
            % cover the range of delays expected from the channel modeling
            %         txWaveform = [txWaveform; zeros(25, 1)]; %#ok
            
            % The initialization time for channel modeling is set each subframe
            % to simulate a continuously varying channel
            chcfg.InitTime = subframeNo/1000;
            
            % Pass data through channel model
            %         rxWaveform = lteFadingChannel(chcfg, txWaveform);
            
            % Add noise at the receiver
            %         v = N*complex(randn(size(rxWaveform)), randn(size(rxWaveform)));
            %         rxWaveform = rxWaveform+v;
            
            
            chcfg.Seed = ueChannelSeed(user);
            
            if (user==1)
                rxWaveform = lteFadingChannel(chcfg,[txWaveform*10^(usersPUSCHpower(user)/20); zeros(25,frc.NTxAnts)]);
                rxWaveform0 = rxWaveform;
            else
                rxWaveform = rxWaveform + ...
                    lteFadingChannel(chcfg,[txWaveform*10^(usersPUSCHpower(user)/20); zeros(25,frc.NTxAnts)]);
            end;
            
        end; % end of user
        
        % Add Noise at the receiver
        noise = N*complex(randn(size(rxWaveform)),randn(size(rxWaveform)));
        rxWaveform = rxWaveform + noise;
        rxWaveform = rxWaveform*sqrt(double(info.Nfft));    % power normalize
        
%         figure();
%         plot((1:size(rxWaveform,1)),abs(rxWaveform0(:,1)),'r');grid on;
%         figure();
%         plot((1:size(rxWaveform,1)),abs(rxWaveform(:,1)),'b');
%         grid on;
        
        
        
        detindex = 1;
          frc.NCellID = usersNCellID(detindex);
            frc.CyclicShift = userCyclicShift(detindex);
            frc.RNTI = usersRNTI(detindex);
            frc.SeqGroup = userSeqGroup(detindex); 
            
          % load dsp data,unit test
        if 1  
        offsetused = 0;
%         dspdata = load('E:\刘新伟拷入文件\cell_search_40subframe_Ant0_Iq_data子帧7.txt');
%         rxWaveform = [[dspdata(:,1) + 1j*dspdata(:,2)];zeros(25,1)];
%         dspdata = load('E:\刘新伟拷入文件\20150114\SF7_TD_10RB_n50.mat');
%         rxWaveform = [dspdata.Msf;zeros(25,1)];
        rxWaveform0 = rxWaveform;
        
        drs_pos = 2048*3+160+144*3+1;
        drs_index = [drs_pos:drs_pos+2048-1 drs_pos+15360:drs_pos+2048-1+15360];
        drs_td = rxWaveform0(drs_index);
        P_td = mean(drs_td.*conj(drs_td));    
        P_td_dB = 10*log10(P_td);    
        
        figure();
        plot(10*log10(drs_td.*conj(drs_td)));
        grid on;
        hold on;
        plot([1:4096],P_td_dB,'r-');
         
        end;
        
        
        % Calculate synchronization offset
        offset = lteULFrameOffset(frc, frc.PUSCH, rxWaveform);    % 利用综测仪信号验证，这个估计还是对星座点有改进的
        if (offset < 25)
            offsetused = offset;
        end;
                     
           offsetused = 24;
%         load('E:\刘新伟拷入文件\FD_data.mat');
%         rxSubframe222 = lteSCFDMADemodulate(frc, ...
%             rxWaveform(1+offsetused:end, :));
%         subdata = rxSubframe222 - rxSubframe0;
%         plot(abs(subdata));
        
        % SC-FDMA demodulation
        rxSubframe = lteSCFDMADemodulate(frc, ...
            rxWaveform(1+offsetused:end, :));
         rxSubframe0 = lteSCFDMADemodulate(frc, ...
            rxWaveform0(1+offsetused:end, :));
       
        % power normalize
        rxSubframe = rxSubframe./sqrt(2048);
        rxSubframe0 = rxSubframe0./sqrt(2048);
        figure();
        plot(abs(rxSubframe0(:,4)));
        grid on;
        
        %
        
        
        % Channel and noise power spectral density estimation
        [estChannelGrid, noiseEst] = lteULChannelEstimate(frc, ... 
            frc.PUSCH, cec, rxSubframe);
         [estChannelGrid0, noiseEst0] = lteULChannelEstimate(frc, ... 
            frc.PUSCH, cec, rxSubframe0);
       
        % Extract resource elements (REs) corresponding to the PUSCH from
        % the given subframe across all receive antennas and channel
        % estimates
        puschIndices = ltePUSCHIndices(frc, frc.PUSCH);
        [puschRx, puschEstCh] = lteExtractResources( ...
            puschIndices, rxSubframe, estChannelGrid);
         [puschRx0, puschEstCh0] = lteExtractResources( ...
            puschIndices, rxSubframe0, estChannelGrid0);
        
        % bj: calc power               
        H = estChannelGrid(121:120+12*1,[4 11]);
        H0 = estChannelGrid0(121:120+12*1,[4 11]);
%         figure();
%         plot(real(H(:,1)));
         Sp(subframeNo+1,1) = mean(mean(H.*conj(H)))/chcfg.NRxAnts;
        Sp0(subframeNo+1,1) = mean(mean(H0.*conj(H0)))/chcfg.NRxAnts;
        Sp_TD_liner = Sp*12/2048;
        Sp_TD_dB = 10*log10(Sp_TD_liner);
        Sp_TD_BB_dB = Sp_TD_dB+33;
                
        % MMSE equalization
        rxSymbols = lteEqualizeMMSE(puschRx, puschEstCh, noiseEst);
        
%         figure();
%         scatter(real(rxSymbols),imag(rxSymbols));

        % Update frc.PUSCH to carry complete information of the UL-SCH
        % coding configuration
        % bj: dsp data decode, unit test
%        
%         harqProc(harqIdx).trBlkSize = 4584;
        
        frc.PUSCH = lteULSCHInfo(frc, ...
            frc.PUSCH, harqProc(harqIdx).trBlkSize, 'chsconcat');

        % Decode the PUSCH        
        rxEncodedBits = ltePUSCHDecode(frc, frc.PUSCH, rxSymbols);
               
        % Decode the UL-SCH channel and store the block CRC error for given
        % HARQ process harqIdx
        trBlkSize = trBlkSizes(mod(subframeNo, 10)+1);
        
        % bj:unit test
%         trBlkSize =  4584;
%          frc.PUSCH.RVSeq = 1;
%         frc.PUSCH.RV = 0;
        
        [rxDecodedBits, harqProc(harqIdx).crc, ...
            harqProc(harqIdx).decState] = lteULSCHDecode(...
            frc, frc.PUSCH, trBlkSize, ...
            rxEncodedBits, harqProc(harqIdx).decState);

        % Store the CRC calculation and total number of bits per subframe
        % successfully decoded
        blkCRC(subframeNo+1) = harqProc(harqIdx).crc;
        bitTp(subframeNo+1) = ...
            harqProc(harqIdx).trBlkSize.*(1-harqProc(harqIdx).crc);
        

    end;   % end of subframe     
    
    % Record the block CRC error and bit throughput for the total number of
    % frames simulated at a particular SNR
    totalBLKCRC(resultIndex, :) = blkCRC;
    bitThroughput(resultIndex, :) = bitTp;
    resultIndex = resultIndex + 1;
    
    % stat Sp and Sp0
    Spave(nSNR) = mean(Sp);
    Spstd(nSNR) = std(Sp,1,1);
    Sp0ave(nSNR) = mean(Sp0);
    Sp0std(nSNR) = std(Sp0,1,1);
    
    
end; % end of snrdB;


if 0
    figure();
    plot(SNRIn,Spave,'b-o',SNRIn,Sp0ave,'r-*');
    grid on;
    legend('real','ideal');
    xlabel('SNR (dB)');
    ylabel('liner mean power');
    title('PUSCH received average power');
    
    figure();
    plot(SNRIn,Spstd,'b-o',SNRIn,Sp0std,'r-*');
    grid on;
    legend('real','ideal');
    xlabel('SNR (dB)');
    ylabel('liner Std');
    title('PUSCH received Std');
    
end;



%% Display Throughput Results
% The throughput results are plotted as a percentage of total capacity and
% actual bit throughput for the range of SNR values input using
% <matlab:edit('hPUSCHResults.m') hPUSCHResults.m>.

% Throughput calculation as a percentage
throughput = 100*(1-mean(totalBLKCRC, 2)).';

residualBler = 1 - throughput/100;
figure();
semilogy(SNRIn,residualBler,'r-o');
 title([ num2str(length(rvSeq))  'Tx residualBler for ', num2str(numSubframes) ' SubFrames']);
    xlabel('SNR (dB)'); ylabel('residualBler');
grid on;
set(gca,  'XLim', [min(SNRIn) max(SNRIn)],'YLim', [0.001 1]);

hPUSCHResults(SNRIn, numSubframes, trBlkSizes, throughput, bitThroughput);

%% Appendix
% This example uses the helper functions:
%
% * <matlab:edit('hPUSCHHARQScheduling.m') hPUSCHHARQScheduling.m>
% * <matlab:edit('hPUSCHNewHARQProcess.m') hPUSCHNewHARQProcess.m>
% * <matlab:edit('hPUSCHResults.m') hPUSCHResults.m>

%% Selected Bibliography
% # 3GPP TS 36.104

displayEndOfDemoMessage(mfilename)
