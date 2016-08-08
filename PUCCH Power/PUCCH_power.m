%% PUCCH1a Multi User ACK Missed Detection Probability Conformance Test
% This example shows how to use the LTE System Toolbox(TM) to measure the
% probability of Acknowledgment (ACK) missed detection for multiuser
% Physical Uplink Control Channel (PUCCH) 1a. The test conditions are
% defined in TS36.104 Section 8.3.4.1 [ <#8 1> ].

% Copyright 2011-2014 The MathWorks, Inc.

%% Introduction
% In this example, four different UEs are configured, each of which
% transmits a PUCCH format 1a signal. Appropriate Demodulation Reference
% Signals (DRS) are also generated. For each considered SNR value, the
% transmitted signals are fed through different channels and added,
% together with Gaussian noise. This simulates the reception of the signals
% from four different UEs at a base station. The receiver decodes the
% PUCCH1a for the user of interest and the probability of ACK missed
% detection is measured. This example uses a simulation length of 10
% subframes. This value has been chosen to speed up the simulation. A
% larger value should be chosen to obtain more accurate results. The target
% defined in TS36.104 Section 8.3.4.1 [ <#8 1> ] for 1.4 MHz bandwidth (6
% Resource Blocks-RBs) and a single transmit antenna is an ACK missed
% detection probability not exceeding 1% at an SNR of -4.1 dB. The test is
% defined for 1 transmit antenna.
clear;

numSubframes = 10;                          % Number of subframes
SNRdB = [-10.1 -8.1 -6.1 -4.1 -2.1];        % SNR range
% SNRdB = [-20 -15 -10];    
NTxAnts = 1;                                % Number of transmit antennas

%% UE 1 Configuration
% Create a User Equipment (UE) configuration structure. These parameters
% are common for all the users.


ue = struct;                  % UE config structure
ue.NULRB = 100;                 % 6 resource blocks (1.4 MHz)
ue.CyclicPrefixUL = 'Normal'; % Normal cyclic prefix
ue.Hopping = 'Off';         % No frequency hopping
% ue.NCellID = 150;           % Cell id as specified in TS36.104 Appendix A9
ue.Shortened = 0;           % No SRS transmission
ue.NTxAnts = NTxAnts;

%% PUCCH 1a Configuration
% We intend to transmit an ACK via a PUCCH of Format 1, so we create an
% appropriate configuration structure |pucch|. We give the cell an 
% arbitrary identification number, and set up PUCCH Resource Indices, 
% transmission powers and channel seeds for each user. A different 
% random channel seed for each user ensures that each experiences 
% different channel conditions. 

% Hybrid Automatic Repeat Request (HARQ)  indicator bit set to one. Only
% one bit is required for PUCCH 1a
ACK = 1;         
pucch = struct;  % PUCCH config structure
% Set the size of resources allocated to PUCCH format 2. This affects the
% location of PUCCH 1 transmission
pucch.ResourceSize = 0; 
% Delta shift PUCCH parameter as specified in TS36.104 Appendix A9 [ <#8 1> ]
pucch.DeltaShift = 2;                       
% Number of cyclic shifts used for PUCCH format 1 in resource blocks with a
% mixture of formats 1 and 2. This is the N1cs parameter as specified in
% TS36.104 Appendix A9
pucch.CyclicShifts = 0;   
% Vector of PUCCH resource indices for all UEs as specified in TS36.104
% Appendix A9
% usersPUCCHindices = [2 1 7 14];    
usersPUCCHindices = [2 1 7 14]; 
% PUCCH power for all UEs as specified in TS36.104 Appendix A9
usersPUCCHpower = [0 0 -3 3];
% usersPUCCHpower = [0 6 -3 3];
usersNCellID = [150 150 150 150] ;


%% Propagation Channel Configuration
% This section of the code configures the propagation channels for the four
% UEs. The parameters are specified in the tests described in TS36.104
% Section 8.3.4.1 [ <#8 1> ] and are: ETU 70Hz and 2 receive, i.e. base
% station, antennas is configured to be 2. Each UE will see a different
% channel, therefore a different seed is used in each case. This is
% specified in the |ueChannelSeed| parameter.

channel = struct;                   % Channel config structure
channel.NRxAnts = 2;                % Number of receive antennas
channel.DelayProfile = 'ETU';       % Channel delay profile
channel.DopplerFreq = 70.0;         % Doppler frequency in Hz
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
% The channel estimator is configured using a structure. Here cubic 
% interpolation will be used with an averaging window of 9x9 resource 
% elements. 

cec = struct;        % Channel estimation config structure
cec.TimeWindow = 9;  % Time averaging window size in resource elements
cec.FreqWindow = 9;  % Frequency averaging window size in resource elements
cec.InterpType = 'cubic';         % Cubic interpolation
cec.PilotAverage = 'UserDefined'; % Type of pilot averaging  

%% Simulation Loop for Configured SNR Points
% A loop is used to run the simulation for a set of SNR points, given by
% the vector |SNRdB|. The SNR vector configured here is a range of SNR 
% points including an SNR point at -4.1dB, the SNR at which the test 
% requirement for ACK detection rate (99%) is to be achieved.

% Preallocate memory for missed detection probability vector
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
    
    % Subframe and user loops
    % We now enter two further loops to process multiple subframes and 
    % create each of the users' transmissions. The fading process time 
    % offset, InitTime, is also generated for the current subframe
        
    offsetused = 0;
    for nsf = 1:numSubframes

        % Channel state information: set the init time to the correct value
        % to guarantee continuity of the fading waveform
        channel.InitTime = (nsf-1)/1000;
            
        % Loop for each user   
        for user = 1:4

            % Create resource grid
            ue.NSubframe = mod(nsf-1,10);            
            txgrid = lteULResourceGrid(ue);
            
            % Configure resource index for this user
            pucch.ResourceIdx = usersPUCCHindices(user);
            ue.NCellID  = usersNCellID(user);

            % ACK bit to transmit for the 1st (target) user, the PUCCH
            % Format 1 carries the Hybrid ARQ (HARQ) indicator ACK and for
            % other users it carries a random HARQ indicator. As there is a
            % single indicator, the transmissions will be of Format 1a. The
            % PUCCH Format 1 DRS carries no data.
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
        % The uplink frame timing estimate for UE1 is calculated using
        % the PUCCH 1 DRS signals and then used to demodulate the             
        % SC-FDMA signal.
        % An offset within the range of delays expected from the channel 
        % modeling (a combination of implementation delay and channel 
        % delay spread) indicates success.
        offset = lteULFrameOffsetPUCCH1(ue,pucch,rxwave);
        if (offset<25)
            offsetused = offset;
        end
        
        % SC-FDMA demodulation
        % The resulting grid (rxgrid) is a 3-dimensional matrix. The number
        % of rows represents the number of subcarriers. The number of
        % columns equals the number of SC-FDMA symbols in a subframe. The
        % number of subcarriers and symbols is the same for the returned
        % grid from lteSCFDMADemodulate as the grid passed into
        % lteSCFDMAModulate. The number of planes (3rd dimension) in the
        % grid corresponds to the number of receive antenna.
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

figure();
plot(SNRdB,Spave,'b-o',SNRdB,Sp0ave,'r-*');
grid on;
legend('real','ideal');
xlabel('SNR (dB)');
ylabel('liner mean power');
title('PUCCH 1a received average power');  

figure();
plot(SNRdB,Spstd,'b-o',SNRdB,Sp0std,'r-*');
grid on;
legend('real','ideal');
xlabel('SNR (dB)');
ylabel('liner Std');
title('PUCCH 1a received Std');  



%% Results
% Finally we plot the simulated results against the target performance as 
% stipulated in the standard.
figure();
semilogy(SNRdB,PMISS,'b-o','LineWidth',2,'MarkerSize',7); 
grid on;
hold on;
semilogy(-4.4,0.01,'rx','LineWidth',2,'MarkerSize',7);
xlabel('SNR (dB)');
ylabel('Probability of ACK missed detection');
title('Multi user PUCCH Format 1a test (TS36.104 Section 8.3.4.1)');   
axis([SNRdB(1)-0.1 SNRdB(end)+0.1 0.001 1.0]);
legend('simulated performance','target');

%% Selected Bibliography
% # 3GPP TS 36.104.

displayEndOfDemoMessage(mfilename);


close all;
