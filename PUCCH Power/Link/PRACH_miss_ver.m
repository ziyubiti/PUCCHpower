%% PRACH Detection 

clear;
close all;
clc;
path(path,'../Model');

numSubframes = 100000;  % Number of subframes frames to simulate at each SNR
SNRdB = [-15.0:1.0:0.0];  % SNR points to simulate
foffset = 270.0;                        % Frequency offset in Hertz

%% UE Configuration
% User Equipment (UE) settings are specified in the structure |ue|.

ue.NULRB = 100;                   % 6 Resource Blocks
ue.DuplexMode = 'TDD';          % Frequency Division Duplexing (FDD)
ue.CyclicPrefixUL = 'Normal';   % Normal cyclic prefix length
ue.NTxAnts = 1;                 % Number of transmission antennas
ue.TDDConfig = 2;
ue.SSC = 7;
ue.CyclicPrefix = 'Normal';
 
    
%% PRACH Configuration
prach.Format = 4;          % PRACH format: TS36.104, Table 8.4.2.1-1
prach.SeqIdx = 0;%22;         % Logical sequence index: TS36.141, Table A.6-1
prach.CyclicShiftIdx = 4; %1;  % Cyclic shift index: TS36.141, Table A.6-1
prach.HighSpeed = 0;       % Normal mode: TS36.104, Table 8.4.2.1-1
prach.FreqOffset = 0;      % Default frequency location
prach.PreambleIdx = 0; %32;    % Preamble index: TS36.141, Table A.6-1
prach.ConfigIdx = 51;
prach.FreqIdx = 0;

info = ltePRACHInfo(ue, prach);  % PRACH information
    
%% Propagation Channel Configuration

chcfg.NRxAnts = 2;                       % Number of receive antenna
chcfg.DelayProfile = 'ETU'; %'ETU';              % Delay profile
chcfg.DopplerFreq = 70.0;                % Doppler frequency  
chcfg.MIMOCorrelation = 'Low';           % MIMO correlation
chcfg.Seed = 1;                          % Channel seed   
chcfg.NTerms = 16;                       % Oscillators used in fading model
chcfg.ModelType = 'GMEDS';               % Rayleigh fading model type 
chcfg.InitPhase = 'Random';              % Random initial phases        
chcfg.NormalizePathGains = 'On';         % Normalize delay profile power
chcfg.NormalizeTxAnts = 'On';            % Normalize for transmit antennas
chcfg.SamplingRate = info.SamplingRate;  % Sampling rate
    
%% Loop for SNR Values

% Initialize the random number generator stream
rng('default');
    
% Initialize variables storing probability of detection at each SNR
pDetection = zeros(size(SNRdB));

for nSNR = 1:length(SNRdB)
    fprintf('\nSimulating at %g dB SNR for a total %d SubFrame(s)', ...
        SNRdB(nSNR), numSubframes);
    
    % Scale noise to ensure the desired SNR after SC-FDMA demodulation
    ulinfo = lteSCFDMAInfo(ue);
    SNR = 10^(SNRdB(nSNR)/20);        
    N = 1/(SNR*sqrt(double(ulinfo.Nfft)))/sqrt(2.0); 

    % Detected preamble count
    detectedCount = 0;  
    
    % Loop for each subframe
    for nsf = 1:numSubframes

        % PRACH transmission
        ue.NSubframe = mod(nsf-1, 10);
        ue.NFrame = fix((nsf-1)/10);
        
        if (ue.NSubframe ~= 1)
            continue;
        end;
        
        % Set PRACH timing offset in us as per TS36.141, Figure 8.4.1.4.2-2
        prach.TimingOffset = info.BaseOffset + ue.NSubframe/10.0;
       
        % Generate transmit wave
        txwave = ltePRACH(ue, prach);             

        % Channel modeling
        chcfg.InitTime = (nsf-1)/1000;
        [rxwave, fadinginfo] = lteFadingChannel(chcfg, ...
                                [txwave; zeros(25, 1)]);

        % Add noise
        noise = N*complex(randn(size(rxwave)), randn(size(rxwave)));            
        rxwave = rxwave + noise;            

        % Remove the implementation delay of the channel modeling
        rxwave = rxwave((fadinginfo.ChannelFilterDelay + 1):end, :);  

        % Apply frequency offset
        t = ((0:size(rxwave, 1)-1)/chcfg.SamplingRate).';
        rxwave = rxwave .* repmat(exp(1i*2*pi*foffset*t), ...
            1, size(rxwave, 2));

        % PRACH detection for all cell preamble indices
        [detected, offsets] = ltePRACHDetect(ue, prach, rxwave, (0:63).');

        % Test for preamble detection
        if (length(detected)==1)
            
            % Test for correct preamble detection
            if (detected==prach.PreambleIdx)         
                
                % Calculate timing estimation error. The true offset is
                % PRACH offset plus channel delay
                trueOffset = prach.TimingOffset/1e6 + 310e-9;
                measuredOffset = offsets(1)/chcfg.SamplingRate;
                timingerror = abs(measuredOffset-trueOffset);
                
                % Test for acceptable timing error
                if (timingerror<=2.08e-6)
                    detectedCount = detectedCount + 1; % Detected preamble
                else
                    disp('Timing error');
                end
            else
%                 disp('Detected incorrect preamble');
                ;
            end
        else
%             disp('Detected multiple or zero preambles');
            ;
        end

    end % of subframe loop

    % Compute final detection probability for this SNR
    pDetection(nSNR) = detectedCount/numSubframes;

end % of SNR loop

%% Analysis

% hPRACHDetectionResultsV2(SNRdB, numSubframes, pDetection);    % for FDD;
hPRACHDetectionResultsV2(SNRdB, numSubframes/10, pDetection*10);   % for TDD      


% busrt 0,record etu70
snr = [-20,-19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7];
pd = [0.214800000000000,0.316500000000000,0.431700000000000,0.552000000000000,0.661600000000000,...
    0.760000000000000,0.834100000000000,0.898800000000000,0.941000000000000,0.966300000000000,0.983100000000000,...
    0.992200000000000,0.996700000000000,0.999100000000000];
% awgn
snr = [-24,-23,-22,-21,-20,-19,-18,-17,-16,-15];
pd = [0.00830000000000000,0.0283000000000000,0.0781000000000000,0.181000000000000,0.377700000000000,0.634300000000000,0.866200000000000,0.974800000000000,0.997900000000000,1];


% burst 4,etu70

% awgn
snr = [-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6];
pd = [0.000740000000000000,0.00208000000000000,0.00585000000000000,...
    0.0156900000000000,0.0329400000000000,0.0593600000000000,0.0847700000000000,0.0972300000000000,0.0997800000000000,0.100000000000000,0.100000000000000];



SNRdB = snr;
pDetection = pd;


displayEndOfDemoMessage(mfilename)