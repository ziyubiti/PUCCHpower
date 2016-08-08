%% PUCCH1a Multi User ACK Missed Detection Probability 

clear;
close all;
clc;

numSubframes = 10000;                          % Number of subframes
% SNRdB = [-10.1 -8.1 -6.1 -4.1 -2.1];        % SNR range
SNRdB = [ -10.0:1.0:-6.0 ];    
NTxAnts = 1;                                % Number of transmit antennas

%% UE 1 Configuration,common for all the users.
ue = struct;                  % UE config structure
ue.NULRB = 100;                 % 
ue.CyclicPrefixUL = 'Normal'; % Normal cyclic prefix
ue.Hopping = 'Off';         % No frequency hopping
% ue.NCellID = 150;           % Cell id 
ue.Shortened = 0;           % No SRS transmission
ue.NTxAnts = NTxAnts;

%% PUCCH 1a Configuration
ACK = 1;         
pucch = struct;  % PUCCH config structure
% Set the size of resources allocated to PUCCH format 2. This affects the
% location of PUCCH 1 transmission
pucch.ResourceSize = 0; 
pucch.DeltaShift = 2;                       
pucch.CyclicShifts = 0;   
% usersPUCCHindices = [2 1 7 14];    
usersPUCCHindices = [2 1 7 14]; 
usersPUCCHpower = [0 0 -3 3];
% usersPUCCHpower = [0 6 -3 3];
usersNCellID = [150 150 150 150] ;
NumUser = 1;             % user number


%% Propagation Channel Configuration
channel = struct;                   % Channel config structure
channel.NRxAnts = 2;                % Number of receive antennas
channel.DelayProfile = 'EPA';       % Channel delay profile
channel.DopplerFreq = 5.0;         % Doppler frequency in Hz
channel.MIMOCorrelation = 'Low';    % Low MIMO correlation
channel.NTerms = 16;                % Oscillators used in fading model
channel.ModelType = 'GMEDS';        % Rayleigh fading model type    
channel.InitPhase = 'Random';       % Random initial phases     
channel.NormalizePathGains = 'On';  % Normalize delay profile power   
channel.NormalizeTxAnts = 'On';     % Normalize for transmit antennas

% SC-FDMA modulation information: required to get the sampling rate
info = lteSCFDMAInfo(ue);
channel.SamplingRate = info.SamplingRate;   % Channel sampling rate 

% Channel seeds for each of the 4 UEs (arbitrary)
ueChannelSeed = [1000 100 70 14]; 

%% Channel Estimator Configuration
cec = struct;        % Channel estimation config structure
cec.TimeWindow = 9;  % Time averaging window size in resource elements
cec.FreqWindow = 9;  % Frequency averaging window size in resource elements
cec.InterpType = 'cubic';         % Cubic interpolation
cec.PilotAverage = 'UserDefined'; % Type of pilot averaging  

%% Simulation Loop for Configured SNR Points
PMISS = zeros(size(SNRdB));
Sp = zeros(numSubframes,1);
Sp0 = zeros(numSubframes,1);
Spave = zeros(size(SNRdB));
Spstd = zeros(size(SNRdB));
Sp0ave = zeros(size(SNRdB));
Sp0std = zeros(size(SNRdB));

for nSNR = 1:length(SNRdB)
    
    % Detection failures counter
    missCount = 0;
    falseCount = 0;
    
    % Noise configuration
    SNR = 10^(SNRdB(nSNR)/20);              % Convert dB to linear
    N = 1/(SNR*sqrt(double(info.Nfft)))/sqrt(2.0*ue.NTxAnts);
    rng('default');
    fprintf('\nSimulating at %g dB SNR for a total %d SubFrame(s)', SNRdB(nSNR), numSubframes);
        
    offsetused = 0;
    for nsf = 1:numSubframes

        % Channel state information: set the init time to the correct value
        % to guarantee continuity of the fading waveform
        channel.InitTime = (nsf-1)/1000;
            
        % Loop for each user   
        for user = 1:NumUser

            % Create resource grid
            ue.NSubframe = mod(nsf-1,10);            
            txgrid = lteULResourceGrid(ue);
            
            % Configure resource index for this user
            pucch.ResourceIdx = usersPUCCHindices(user);
            ue.NCellID  = usersNCellID(user);
           
            if (user==1)
                txACK = ACK;
            else
                txACK = randi([0 1],1,1);
            end

            % Generate PUCCH 1 and its DRS        
            % Different users have different relative powers
            pucch1Sym = ltePUCCH1(ue,pucch,txACK)* ...
                10^(usersPUCCHpower(user)/20);    
            pucch1DRSSym = ltePUCCH1DRS(ue,pucch)* ...
                10^(usersPUCCHpower(user)/20); 
            
            % Generate indices for PUCCH 1 and its DRS
            pucch1Indices = ltePUCCH1Indices(ue,pucch);
            pucch1DRSIndices = ltePUCCH1DRSIndices(ue,pucch);
            
            % Map PUCCH 1 and PUCCH 1 DRS to the resource grid
            if (~isempty(txACK))
                txgrid(pucch1Indices) = pucch1Sym;
                txgrid(pucch1DRSIndices) = pucch1DRSSym;
            end
            
            % SC-FDMA modulation
            txwave = lteSCFDMAModulate(ue,txgrid);
            
            % Channel modeling and superposition of received signals.
            % The additional 25 samples added to the end of the waveform
            % are to cover the range of delays expected from the channel
            % modeling (a combination of implementation delay and channel
            % delay spread). On each iteration of the loop we accumulate
            % the sum of each transmitted signal, simulating the reception
            % of all four users at the base station.
            channel.Seed = ueChannelSeed(user);
            if (user==1)
                rxwave = lteFadingChannel(channel,[txwave; zeros(25,NTxAnts)]);
                rxwave0 = rxwave;
            else
                rxwave = rxwave + ...
                    lteFadingChannel(channel,[txwave; zeros(25,NTxAnts)]);
            end;
                        
        end  % end of user
        
         % Add Noise at the receiver
            noise = N*complex(randn(size(rxwave)),randn(size(rxwave)));
            rxwave = rxwave + noise;
            
            
%         figure();
%         plot([1:size(rxwave0)],abs(rxwave0),'r');grid on;
% %         hold on;
%         figure();
%         plot([1:size(rxwave)],abs(rxwave),'b');
%         grid on;
        
        
        
        % Receiver
        
        % Use the resource indices for the user of interest
        detindex = 1;
        pucch.ResourceIdx = usersPUCCHindices(detindex);
        ue.NCellID  = usersNCellID(detindex);
        
        % Synchronization
        offset = lteULFrameOffsetPUCCH1(ue,pucch,rxwave);
        if (offset<25)
            offsetused = offset;
        end
        
        % SC-FDMA demodulation
        rxgrid = lteSCFDMADemodulate(ue,rxwave(1+offsetused:end,:));
        rxgrid0 = lteSCFDMADemodulate(ue,rxwave0(1+offsetused:end,:));
    
        
        % Channel estimation            
        [H,n0] = lteULChannelEstimatePUCCH1(ue,pucch,cec,rxgrid);
        [H0,n00] = lteULChannelEstimatePUCCH1(ue,pucch,cec,rxgrid0);

        % PUCCH 1 indices for UE of interest
        pucch1Indices = ltePUCCH1Indices(ue,pucch);

        % Extract resource elements (REs) corresponding to the PUCCH 1 from
        % the given subframe across all receive antennas and channel
        % estimates
        [pucch1Rx,pucch1H] = lteExtractResources(pucch1Indices,rxgrid,H);
        [pucch1Rx0,pucch1H0] = lteExtractResources(pucch1Indices,rxgrid0,H0);
        
        % bj:calc pucch signal received power
        Sp(nsf,1) = mean(diag(pucch1H*pucch1H'))/channel.NRxAnts;
        Sp0(nsf,1) = mean(diag(pucch1H0*pucch1H0'))/channel.NRxAnts;

        % Minimum Mean Squared Error (MMSE) Equalization
        eqgrid = lteULResourceGrid(ue);    
        eqgrid(pucch1Indices) = lteEqualizeMMSE(pucch1Rx,pucch1H,n0);

        % PUCCH 1 decoding
        rxACK = ltePUCCH1Decode(ue,pucch,1,eqgrid(pucch1Indices));

        % Detect missed (empty rxACK) or incorrect HARQ-ACK (compare
        % against transmitted ACK.
        if (isempty(rxACK) || any(rxACK~=ACK))
            missCount = missCount + 1;  
        end;


    end;   % end nsf
        
    PMISS(nSNR) = missCount/numSubframes;
    
    % stat Sp and Sp0
    Spave(nSNR) = mean(Sp);
    Spstd(nSNR) = std(Sp,1,1);
    Sp0ave(nSNR) = mean(Sp0);
    Sp0std(nSNR) = std(Sp0,1,1);
    
    
end;   % end nSNR

% received power
% figure();
% plot(SNRdB,Spave,'b-o',SNRdB,Sp0ave,'r-*');
% grid on;
% legend('real','ideal');
% xlabel('SNR (dB)');
% ylabel('liner mean power');
% title('PUCCH 1a received average power');  
% 
% figure();
% plot(SNRdB,Spstd,'b-o',SNRdB,Sp0std,'r-*');
% grid on;
% legend('real','ideal');
% xlabel('SNR (dB)');
% ylabel('liner Std');
% title('PUCCH 1a received Std');  



%% Results

figure();
semilogy(SNRdB,PMISS,'b-o','LineWidth',2,'MarkerSize',7); 
grid on;
hold on;
% semilogy(-4.4,0.01,'rx','LineWidth',2,'MarkerSize',7);
xlabel('SNR (dB)');
ylabel('Probability of ACK missed detection');
title('PUCCH Format 1a performance for 10000 subframes ');   
axis([SNRdB(1)-0.1 SNRdB(end)+0.1 0.0001 1.0]);
legend('PUCCH format 1a','target');

displayEndOfDemoMessage(mfilename);

% epa5 record
snr = [-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6];
pm = [0.342000000000000,0.282000000000000,0.219100000000000,0.167900000000000,0.121200000000000,...
    0.0845000000000000,0.0552000000000000,0.0356000000000000,0.0205000000000000,0.0108000000000000,0.00470000000000000,0.00270000000000000,0.00140000000000000];

SNRdB = snr;
PMISS = pm;
