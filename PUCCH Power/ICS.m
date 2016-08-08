%% Cell Search, MIB and SIB1 Recovery 

clear;
close all;
path(path,'Model');
% load eNodeBOutput.mat           % Load I/Q capture of eNodeB output
master_ant0_time_domain_10ms = load('E:\technote20150415\beidongshi\IQ数据3\master_ant1_time_domain_10ms');
eNodeBOutput= master_ant0_time_domain_10ms(1:2:end,1) + 1j*master_ant0_time_domain_10ms(2:2:end,1);

eNodeBOutput = double(eNodeBOutput)/32768; % Scale samples

figure();
plot(abs(eNodeBOutput));
grid on;

% eNodeBOutput = [D2(30721:end,1);D2(1:30720,1)];        
sr = 30.72e6;                   % Sampling rate for loaded samples  

% Set eNodeB basic parameters
enb = struct;                   % eNodeB config structure
enb.DuplexMode = 'FDD';
enb.CyclicPrefix = 'Normal';    
enb.NDLRB = 6;                  % Number of resource blocks
ofdmInfo = lteOFDMInfo(enb);    % Needed to get the sampling rate
downsampling = sr/ofdmInfo.SamplingRate; % Downsampling factor
% Downsample loaded samples
downsampled = resample(eNodeBOutput, 1, downsampling); 
    
%% Cell Search and Synchronization
% Call <matlab:doc('lteCellSearch') lteCellSearch> to obtain the cell 
% identity and timing offset |offset| to the first frame head.
fprintf('\nPerforming cell search...');
[enb.NCellID, offset] = lteCellSearch(enb, downsampled);
% UT
offset = 0;   % UT
enb.NCellID = 136;

downsampled = downsampled(1+offset:end); 
enb.NSubframe = 0;
fprintf('\nCell-wide settings after cell search:\n');
disp(enb);

%% OFDM Demodulation and Channel Estimation  
% The OFDM downsampled I/Q waveform is demodulated to produce a resource
% grid |rgrid|. This is used to perform channel estimation. |hest| is the
% channel estimate, |nest| is an estimate of the noise (for MMSE
% equalization) and |cec| is the channel estimator configuration.
%
% For channel estimation the example assumes 4 cell specific reference
% signals. This means that channel estimates to each receiver antenna from
% all possible cell-specific reference signal ports are available. The true
% number of cell-specific reference signal ports is not yet known. The
% channel estimation is only performed on the first subframe, i.e. using
% the first |L| OFDM symbols in |rxgrid|.

% Channel estimator configuration
cec.PilotAverage = 'UserDefined';     % Type of pilot averaging
cec.FreqWindow = 9;                   % Frequency window size    
cec.TimeWindow = 9;                   % Time window size    
cec.InterpType = 'cubic';             % 2D interpolation type
cec.InterpWindow = 'Centered';        % Interpolation window type
cec.InterpWinSize = 1;                % Interpolation window size  

% Assume 4 cell-specific reference signals for initial decoding attempt;
% ensures channel estimates are available for all cell-specific reference
% signals
enb.CellRefP = 4;   
                    
griddims = lteResourceGridSize(enb); % Resource grid dimensions
L = griddims(2);                     % Number of OFDM symbols in a subframe 
% OFDM demodulate signal 
rxgrid = lteOFDMDemodulate(enb, downsampled);    
% Perform channel estimation
rxgrid(:,[15:end]) = [];    % only sf0
[hest, nest] = lteDLChannelEstimate(enb, cec, rxgrid(:,1:L,:));    
h = (ifft(hest(:,5,1,1)));
figure();
stem(abs(h));
grid on;

sinr_dB = 10*log10(mean(mean(mean(hest.*conj(hest))))/nest);

%% PBCH Demodulation, BCH Decoding, MIB parsing

pbchIndices = ltePBCHIndices(enb);
[pbchRx, pbchHest] = lteExtractResources( ...
    pbchIndices, rxgrid(:,1:L,:), hest(:,1:L,:,:));

% Decode PBCH
[bchBits, pbchSymbols, nfmod4, mib, enb.CellRefP] = ltePBCHDecode( ...
    enb, pbchRx, pbchHest, nest); 

% Parse MIB bits
enb = lteMIB(mib, enb); 

% Incorporate the nfmod4 value output from the function ltePBCHDecode, as
% the NFrame value established from the MIB is the System Frame Number
% (SFN) modulo 4 (it is stored in the MIB as floor(SFN/4))
enb.NFrame = enb.NFrame+nfmod4;

% Display cell wide settings after MIB decoding
fprintf('\nCell-wide settings after MIB decoding:\n');
disp(enb);

%% OFDM Demodulation on Full Bandwidth
% Now that the signal bandwidth is known, OFDM modulation is performed
% using this value.

% Find beginning of frame
offset = lteDLFrameOffset(enb, eNodeBOutput); 
% UT
  offset = 0;
%  offset = 307200*2;  % UT
% Data before the beginning of the frame is not useful
% if first odd frame, extract 2nd frame,offset and FrameNumber are updated
 offset = offset + mod(enb.NFrame,2)*sr/100;
% enb.NFrame = enb.NFrame + 1;
eNodeBOutput = eNodeBOutput(1+offset:end);   

% frame align 10ms

figure();
plot(abs(eNodeBOutput(1:end,1)));
grid on;

% update para,UT;  because 30.72Msps for BW=15M
enb.NDLRB = 100;
enb.CellRefP = 2;
enb.PHICHDuration = 'Normal';
enb.Ng = 'One';
enb.SSC = 7;
enb.TDDConfig = 2;
% OFDM demodulation
rxgrid = lteOFDMDemodulate(enb, eNodeBOutput);   

% load dsp FD IQ data,UT
master_ant0_freq_domain_10ms = load('E:\technote20150415\beidongshi\IQ数据3\master_ant0_freq_domain_10ms');
FD_data= master_ant0_freq_domain_10ms(1:2:end,1) + 1j*master_ant0_freq_domain_10ms(2:2:end,1);
FD_align = reshape(FD_data,[],140);
rxgrid = FD_align([1:1200],:);

figure();
plot(abs(rxgrid(:,14)));
grid on;

% update 15M,UT
% enb.NDLRB = 75;
% enb.PHICHDuration = 'Normal';
% enb.Ng = 'One';
% rxgrid([1:150 1051:1200],:) = []; 

%% SIB1 Decoding
% The following steps are performed in this section
%
% * Physical Control Format Indicator Channel (PCFICH) demodulation, CFI
% decoding
% * PDCCH decoding
% * Blind PDCCH search
% * SIB bits recovery: PDSCH demodulation and DL-SCH decoding
%
% After recovery the SIB CRC should be 0.

% Check this frame contains SIB1
if (mod(enb.NFrame,2)==0)                    

    % Advance to subframe 5 (this is where SIB1 is scheduled) and perform
    % channel estimation
    rxgrid(:,1:(L*5),:) = [];   % Remove subframes 0 to 4
    rxgrid(:,15:end,:) = [];   % dsp unit test, only reserve SF5.
    disp('Performing channel estimation...');
    enb.NSubframe = 5;          % Subframe number is 5
    % Perform channel estimation
    [hest,nest] = lteDLChannelEstimate(enb, cec, rxgrid(:,1:L,:));       
    
    % plot time domain channel pulse
    h = (ifft(hest(:,5,1,1)));
    figure();
    stem(abs(h));
    grid on;
    
    sinr_dB = 10*log10(mean(mean(mean(hest.*conj(hest))))/nest);
    
    % plot freq domain power
    figure();
    plot(abs(reshape(rxgrid,[],1)));
    grid on;

    % PCFICH demodulation, CFI decoding. The CFI is now demodulated and
    % decoded using similar resource extraction and decode functions to
    % those shown already for BCH reception. lteExtractResources is used to
    % extract REs corresponding the PCFICH from the received grid rxgrid
    % and channel estimate hest.
    disp('Decoding CFI...');
    pcfichIndices = ltePCFICHIndices(enb);  % Get PCFICH indices
    [pcfichRx, pcfichHest] = lteExtractResources(pcfichIndices, rxgrid, hest);
    % Decode PCFICH
    cfiBits = ltePCFICHDecode(enb, pcfichRx, pcfichHest, nest);
    enb.CFI = lteCFIDecode(cfiBits);        % Get CFI
    fprintf('Decoded CFI value: %d\n\n', enb.CFI);   
    
    % PDCCH demodulation. The PDCCH is now demodulated and decoded using
    % similar resource extraction and decode functions to those shown
    % already for BCH and CFI reception
    pdcchIndices = ltePDCCHIndices(enb); % Get PDCCH indices
    [pdcchRx, pdcchHest] = lteExtractResources(pdcchIndices, rxgrid, hest);
    % Decode PDCCH and plot constellation, only CHE/EQ/Demod/DeScrambling
    [dciBits, pdcchSymbols] = ltePDCCHDecode(enb, pdcchRx, pdcchHest, nest);
    figure();
    plot(pdcchSymbols,'ko','MarkerFaceColor',[1 0 0], ...
            'MarkerEdgeColor',[0.625 0 0],'MarkerSize',3);
    title('Received PDCCH constellation');    
    grid on;
    

    % PDCCH blind search for System Information (SI) and DCI decoding. The
    % LTE System Toolbox provides full blind search of the PDCCH to find
    % any DCI messages with a specified RNTI, in this case the SI-RNTI.
    disp('PDCCH search for SI-RNTI...');
    pdcch = struct('RNTI', 65535);  
    [dci,dciTb] = ltePDCCHSearchV2(enb, pdcch, dciBits); % Search PDCCH for DCI
    dci = dci{1};
    fprintf('DCI message with SI-RNTI:\n');
    disp(dci);
    % Get the PDSCH configuration from the DCI.
    [pdsch, trblklen] = hPDSCHConfiguration(enb, dci, pdcch.RNTI);
    pdsch.NTurboDecIts = 5;
    fprintf('PDSCH settings after DCI decoding:\n');
    disp(pdsch);

    % PDSCH demodulation and DL-SCH decoding to recover SIB bits.
    % The DCI message is now parsed to give the configuration of the 
    % corresponding PDSCH carrying SIB1, the PDSCH is demodulated and 
    % finally the received bits are DL-SCH decoded to yield the SIB1 bits.

    disp('Decoding SIB1...');        
    % Get PDSCH indices
    pdschIndices = ltePDSCHIndices(enb, pdsch, pdsch.PRBSet);
    [pdschRx, pdschHest] = lteExtractResources(pdschIndices, rxgrid, hest);
    % Decode PDSCH 
    dlschBits = ltePDSCHDecode(enb, pdsch, pdschRx, pdschHest, nest);
    % Decode DL-SCH
    [sib1, crc] = lteDLSCHDecode(enb, pdsch, trblklen, dlschBits);
    fprintf('SIB1 CRC: %d\n',crc);
    
    if crc == 0
        disp('Successful SIB1 recovery.');
    end
    
    sib1_bit = int2str(sib1{1,1}.');
    sib1_bit(find(isspace(sib1_bit))) = [];
%     sib1_bit_hex = dec2hex(bin2dec(sib1_bit));
    
end

%% Appendix
% This example uses the helper function:
%
% * <matlab:edit('hPDSCHConfiguration.m') hPDSCHConfiguration.m>

displayEndOfDemoMessage(mfilename) 
