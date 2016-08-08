%lteULChannelEstimate PUSCH uplink channel estimation
%   [HEST,NOISEEST] = lteULChannelEstimate(...) returns the estimated
%   channel between each transmit and receive antenna and an estimate of
%   the noise power spectral density.
%
%   HEST is an M-by-N-by-NRxAnts-by-NTxAnts array where M is the total
%   number of subcarriers, N is the number of OFDM symbols, NRxAnts is the
%   number of receive antennas and NTxAnts is the number of transmit
%   antennas. Using the reference signals, NOISEEST, provides an estimate
%   of the power spectral density of the noise present on the estimated
%   channel response coefficients. Optionally, the channel estimator can be
%   configured to use the DMRS layers as the reference signal, in which
%   case HEST is an M-by-N-by-NRxAnts-by-NLayers array where NLayers is the
%   number of transmission layers.
%
%   [HEST,NOISEEST] = lteULChannelEstimate(UE,CHS,RXGRID) returns an
%   estimate for the channel by averaging the least squares estimates of
%   the reference symbols across time and copying these across the
%   allocated resource elements within the time frequency grid.
%
%   UE is a structure and must contain the following fields:
%   NULRB          - Number of uplink resource blocks
%   NCellID        - Physical layer cell identity
%   NSubframe      - Subframe number
%   CyclicPrefixUL - Optional. Cyclic prefix length 
%                    ('Normal'(default),'Extended')
%   NTxAnts        - Optional. Number of transmission antennas. 
%                    (1(default),2,4)
%   
%   CHS is a structure and must contain the following fields:
%   PRBSet         - A 1- or 2-column matrix, containing the Physical
%                    Resource Block (PRB) indices corresponding to the
%                    slot wise resource allocations for this PUSCH
%   DynCyclicShift - Optional. Cyclic Shift for DMRS. (0...7) 
%                    (default 0) (yields n2_DMRS)
%   NLayers        - Optional. Number of transmission layers 
%                    (1(default),2,3,4)
%   Only required for NTxAnts=2 or NTxAnts=4:
%      PMI         - Optional. Scalar Precoder Matrix Indication to be used  
%                    during precoding of the DRS reference symbols.
%                    (0(default)...23), depending on NTxAnts and NLayers.
%                    (See <a
%                    href="matlab:help('lteULPMIInfo')">lteULPMIInfo</a>)
%
%   For PRBSet, if a column vector is provided, the resource allocation is
%   the same in both slots of the subframe; the 2-column matrix can be used
%   to specify differing PRBs for each slot in a subframe. Note that the
%   PRB indices are 0-based.
%
%   RXGRID is a 3-dimensional M-by-N-by-NRxAnts array of resource elements.
%   The second dimension of RXGRID contains any number of whole subframes
%   worth of SC-FDMA symbols e.g. for normal cyclic prefix each subframe
%   contains 14 SC-FDMA symbols, therefore N is equal to
%   14*TotalNoSubframes. If TotalNoSubframes>1 then care must be taken to
%   extract the correct subframe from the returned HEST array. The location
%   of the estimated subframe within HEST is specified using the CEC.Window
%   parameter field introduced below.
%
%   [HEST,NOISEEST] = lteULChannelEstimate(UE,CHS,CEC,RXGRID)
%   returns the estimated channel using the method and parameters defined  
%   by the user in the channel estimator configuration structure CEC.
%
%   CEC is a structure which defines the type of channel estimation
%   performed. CEC must contain a set of the following fields:
%   FreqWindow   - Size of window in resource elements used to average
%                  over frequency (odd number, or a multiple of 12)
%   TimeWindow   - Size of window in resource elements used to average
%                  over time (odd number) 
%   InterpType   - Type of 2D interpolation used. See <a href="matlab:
%                  doc('griddata')">griddata</a> for types. 
%                  'None' disables interpolation. 
%   PilotAverage - Optional. Type of pilot averaging. 
%                  ('TestEVM','UserDefined'(default))
%   Reference    - Optional. Specifies point of reference (signals to 
%                  internally generate) for channel estimation. 
%                  ('Antennas'(default), 'Layers', 'None')
%   Only required if RXGRID contains more than one subframe:       
%      Window    - Indicates the position of the subframe in RXGRID for 
%                  which channel estimation is to be performed. Only
%                  channel estimates for this subframe will be calculated,
%                  but estimation will be aided by reference symbols in the
%                  other subframes. The window can be: 'Left', 'Right',
%                  'Centred' or 'Centered'; for 'Centred'/'Centered' the
%                  window size (i.e. the number of subframes) must be odd.
%
%   The 'TestEVM' pilot averaging will ignore other structure fields in
%   CEC, and the method follows that described in TS36.101 Annex F for the
%   purposes of transmitter EVM testing.
%
%   The 'UserDefined' pilot averaging uses a rectangular kernel of size
%   CEC.FreqWindow-by-CEC.TimeWindow and performs a 2D filtering operation
%   upon the pilots. Note that pilots near the edge of the resource grid
%   will be averaged less as they have no neighbors outside of the grid.
%   For CEC.FreqWindow = 12*X (i.e. any multiple of 12) and CEC.TimeWindow
%   = 1 the estimator will enter a special case where an averaging window
%   of (12*X)-in-frequency will be used to average the pilot estimates; the
%   averaging is always applied across (12*X) subcarriers, even at the
%   upper and lower band edges; therefore the first (6*X) symbols at the
%   upper and lower band edge will have the same channel estimate. This
%   operation ensures that averaging is always done on 12 (or a multiple of
%   12) symbols. This provides the appropriate "despreading" operation
%   required for the case of Uplink MIMO where the DMRS signals associated
%   with each layer occupy the same time/frequency locations but use
%   different orthogonal covers to allow them to be differentiated at the
%   receiver.
%
%   Reference='Antennas' will use the PUSCH DMRS after precoding onto the
%   transmission antennas as the reference for channel estimation. In this
%   case, the precoding matrix indicated in CHS.PMI is used to precode the
%   DMRS layers onto antennas, and the channel estimate HEST will be of
%   size M-by-N-by-NRxAnts-by-NTxAnts where NTxAnts is given by
%   CHS.NTxAnts. Reference='Layers' will use the PUSCH DMRS without
%   precoding as the reference for channel estimation. The channel estimate
%   HEST will be of size M-by-N-by-NRxAnts-by-NLayers where NLayers is
%   given by CHS.NLayers. Reference='None' will generate no internal
%   reference signals, and the channel estimation can be performed on
%   arbitrary known REs as given by the REFGRID argument below. This
%   approach can be used to provide a REFGRID containing for example the
%   SRS signals created on all NTxAnts, allowing for full-rank channel
%   estimation for the purposes of PMI selection when the PUSCH is
%   transmitted with less than full rank.
%
%   For CEC.InterpType='None', no interpolation between pilot symbol
%   estimates will be performed and no virtual pilots will be created. HEST
%   will contain channel estimates in the locations of transmitted
%   reference symbols for each received antenna, and all other elements of
%   HEST will be zero. Note that the averaging of pilot symbol estimates
%   described by CEC.TimeWindow and CEC.FreqWindow will still be performed.
%
%   In the case where RXGRID contains more than one subframe, CEC.Window
%   provides control of the location of the subframe for which channel
%   estimation is to be performed. This allows channel estimation for the
%   subframe of interest to be aided by the presence of pilot symbols in
%   the subframes before and/or after that subframe. For example, if RXGRID
%   contains 5 subframes, 'Left' will estimate the first subframe in
%   RXGRID, 'Centred'/'Centered' will estimate the third (middle) subframe
%   and 'Right' will estimate the last subframe. Note that the parameter
%   ENB.NSubframe corresponds to the chosen subframe, so for example with 3
%   subframes and Window='Right', RXGRID corresponds to subframes
%   ENB.NSubframe-2 ... ENB.NSubframe. The HEST output will be the same
%   size as RXGRID and will correspond to the same subframe numbers. All
%   locations other than the estimated subframe will contain zeros.
%
%   [HEST,NOISEEST] = lteULChannelEstimate(UE,CHS,CEC,RXGRID,REFGRID)
%   returns the estimated channel using the method and parameters defined
%   by the channel estimation configuration structure and the additional
%   information about the transmitted symbols found in REFGRID.
%
%   REFGRID is a 3-dimensional M-by-N-by-NTxAnts array containing known
%   transmitted data symbols in their correct locations. All other
%   locations i.e. DRS Symbols and unknown data symbol locations should be
%   represented by a NaN. RXGRID and REFGRID must have the same dimensions.
%   A typical use for REFGRID is to provide values of the SRS transmitted
%   at some point during the time span of RXGRID, in order that the channel
%   estimation can be enhanced. For CEC.InterpType='None', values in
%   REFGRID are treated as reference symbols and the resulting HEST will
%   contain non-zero values in their locations. 
%
%   [HEST,NOISEEST] = lteULChannelEstimate(UE,CHS,RXGRID,REFGRID) 
%   returns the estimated channel using the estimation method as described 
%   in TS36.101 Annex F4. The method described utilizes extra channel
%   information obtained through information of the transmitted symbols
%   found in REFGRID. This additional information allows for an improved
%   estimate of the channel and is required for accurate EVM measurements. 
%   RXGRID and REFGRID must only contain a whole subframe worth of OFDM
%   symbols i.e. for normal cyclic prefix each subframe contains 14 OFDM 
%   symbols, therefore N is equal to 14. 
%
%   Example: 
%   Estimate the channel characteristics given the received resource grid.
%
%   ue = lteRMCUL('A3-2');
%   ue.TotSubframes = 1;
%   cec = struct('FreqWindow',7,'TimeWindow',1,'InterpType','cubic');
%   txWaveform = lteRMCULTool(ue,[1;0;0;1]);
%   rxWaveform = txWaveform;
%   rxGrid = lteSCFDMADemodulate(ue,rxWaveform);
%   hest = lteULChannelEstimate(ue,ue.PUSCH,cec,rxGrid);
%  
%   See also lteSCFDMADemodulate, lteULFrameOffset, lteEqualizeULMIMO,
%   lteEqualizeMMSE, lteEqualizeZF, lteULPerfectChannelEstimate, griddata.

%   Copyright 2010-2014 The MathWorks, Inc.


function [Hest, NoisePowerEst, Hls] = lteULChannelEstimateV2(varargin)

    ue = varargin{1};
    pusch = varargin{2};
    
    if (~isfield(ue,'NTxAnts'))
        ue.NTxAnts = 1;
        defaultValueWarning('NTxAnts','1'); 
    end
    
    if(isstruct(varargin{3}))
        cec = varargin{3};
        if ~isfield(cec,'PilotAverage')
            cec.PilotAverage = 'UserDefined';  
            defaultValueWarning('PilotAverage','''UserDefined''');
        end
        if (~isfield(cec,'Reference'))
            cec.Reference = 'Antennas';
            defaultValueWarning('Reference','''Antennas''');
        end
        if (strcmpi(cec.Reference,'Antennas')==1)
            NTx = ue.NTxAnts;
            refString = 'UE.NTxAnts';
        elseif (strcmpi(cec.Reference,'Layers')==1)
            NTx = pusch.NLayers;
            refString = 'CHS.Nlayers';
        else
            if (strcmpi(cec.Reference,'None')~=1)
                error('lte:error','Reference must be ''Antennas'', ''Layers'' or ''None'', see help for details.');
            end            
        end
        rxGrid = varargin{4};        
        if nargin==5
            refGrid = varargin{5};
            if (strcmpi(cec.Reference,'None')==1)
                NTx = size(refGrid,3);
            end
        else   
            if (strcmpi(cec.Reference,'None')==1)
                % special case, no internal reference signal is enabled
                % and no REFGRID is provided.
                error('lte:error','If Reference=''None'' is used, a reference grid input must be provided.');                
            end
            refGrid = NaN([size(rxGrid,1) size(rxGrid,2) NTx]);
        end
    else
        % If no configuration structure then use TestEvm method
        rxGrid = varargin{3};
        cec.PilotAverage = 'TestEVM';    
        cec.Reference = 'Antennas';
        refString = 'UE.NTxAnts';
        NTx = ue.NTxAnts;
        if nargin==4
            refGrid = varargin{4};
        else
            refGrid = NaN([size(rxGrid,1) size(rxGrid,2) NTx]);
        end
    end   
    if (size(refGrid,3) ~= NTx)
        error('lte:error','The reference grid input must contain the same number of planes (size of the 3rd dimension) as the configured reference signal; ''%s'' is the configured reference signal therefore %s=%d planes are required.',cec.Reference,refString,NTx);
    end
    
    % Get dimensions of resource grid
    Dims = lteULResourceGridSize(ue);
    K = Dims(1);
    Lsf = Dims(2);
   
    % Initialize size of estimated channel grid
    Hls = zeros(3,2400);         % UT,temp added.
    
    Hest = zeros([size(rxGrid,1) size(rxGrid,2) size(rxGrid,3) NTx]);
    
    % Determine number of Rx Antennas
    NRxAnts = size(rxGrid,3);     
    
    % Preallocate noise power estimate vector for speed
    noiseVec = zeros(size(NRxAnts,NTx));        
    
    if (~isempty(rxGrid))
        
        for rxANT = 1:NRxAnts
            for tx = 1:NTx
                % Extract pilot symbols from received grid and calculate least
                % squares estimate
                [ls_estimates,focalsf,refIndices,pilotSym]= GetPilotEstimates(ue,pusch,tx-1,rxGrid(:,:,rxANT),cec,refGrid(:,:,tx));                
                Hls = ls_estimates;
                if strcmpi(cec.PilotAverage,'TestEVM')
                    % Extract pilot symbols from first slot of subframe,
                    % average and apply to rest of PUSCH region
                    vec = ls_estimates(:,ls_estimates(2,:)<Lsf/2);
                    subcarriers = unique(vec(1,:));
                    Xarr = zeros(numel(subcarriers),1);
                    for i = 1:numel(subcarriers)
                        Xarr(i) = mean(vec(3,(vec(1,:)==subcarriers(i))));
                    end                

                    firstsc_s1 = min(vec(1,:));
                    lastsc_s1 = max(vec(1,:));
                    Hest(firstsc_s1:lastsc_s1,1+Lsf*(focalsf-1):Lsf/2+Lsf*(focalsf-1),rxANT,tx) = repmat(Xarr,1,Lsf/2);

                    % Extract pilot symbols from first slot of subframe,
                    % average and apply to rest of PUSCH region
                    vec = ls_estimates(:,ls_estimates(2,:)>Lsf/2);
                    subcarriers = unique(vec(1,:));
                    for i = 1:numel(subcarriers)
                        Xarr(i) = mean(vec(3,(vec(1,:)==subcarriers(i))));
                    end                

                    firstsc_s2 = min(vec(1,:));
                    lastsc_s2 =  max(vec(1,:));
                    Hest(firstsc_s2:lastsc_s2,1+Lsf/2+Lsf*(focalsf-1):Lsf/2*2+Lsf*(focalsf-1),rxANT,tx) = repmat(Xarr,1,Lsf/2);

                elseif ~(strcmpi(cec.PilotAverage,'UserDefined'))
                    error('lte:error','PilotAverage must be "UserDefined" or "TestEVM", see help for details.');

                else

                    % Average ls_estimates to reduce noise
                    [p_est, scalingVec]= PilotAverage(cec,Hest,ls_estimates);

                    % create virtual pilots:
                    if (~strcmpi(cec.InterpType,'None'))
                        if (~isempty(pilotSym))
                            % When DMRS reference is used:

                            % Create virtual pilots for the focal subframe only            
                            rpilots = [p_est(:,p_est(2,:)==pilotSym(1)) p_est(:,p_est(2,:)==pilotSym(2))];

                            % Create virtual pilots to aid interpolation. VPs are
                            % created in a fashion optimized for PUSCH DMRS.
                            vps = CreateVirtualPilot(ue,rpilots,focalsf);

                        else
                            % When only reference grid is used:
                            % Create virtual pilots to aid interpolation. VPs are
                            % created around the grid edge as we do not know the
                            % relationship between where the reference information
                            % is and where the area of subsequent demodulation is. 
                            % With VPs around the grid edge, we can return some sort 
                            % of estimate at any location in the grid. 
                            vps = createEdgeVirtualPilots(ue,p_est);
                        end

                        % Determine and null out VPs if they conflict with any 
                        % estimates provided by reference grid of known data
                        % symbols
                        if size(refGrid,2)/Lsf>1

                            local_dims = size(refGrid);
                            % Find any virtual pilots created outwith the bounds of
                            % the resource grid and store them in a vector
                            outRangeVPs = [vps(:,vps(1,:)<min(p_est(1,:))) vps(:,vps(1,:)>max(p_est(1,:))) vps(:,vps(2,:)<1) vps(:,vps(2,:)>local_dims(2))];                    
                            % Null the pilots lying outside the resource grid
                            % location
                            vps(:,vps(1,:)<min(p_est(1,:))) = [];
                            vps(:,vps(1,:)>max(p_est(1,:))) = [];
                            vps(:,vps(2,:)<1) = []; 
                            vps(:,vps(2,:)>local_dims(2)) = [];
                            % Determine linear indices of virtual pilots found
                            % within resource grid allocation
                            [vpind] = sub2ind([local_dims(1),local_dims(2)],vps(1,:),vps(2,:));
                            % Create grid of NaNs of size refGrid
                            pilotGrid = NaN([local_dims(1),local_dims(2)]);
                            % Place these VPs into dummy grid
                            pilotGrid(vpind) = vps(3,:);
                            % Place known data symbol least squares estimates into
                            % dummy grid locations, overwriting any VPs that were
                            % created in their locations
                            pilotGrid(refIndices) = p_est(3,:);
                            %Extract all desired symbols from the dummy grid
                            [sc,sym] = ind2sub([local_dims(1),local_dims(2)], find(~isnan(pilotGrid)));                    
                            p_est = [double(sc).' ; double(sym).' ; pilotGrid(~isnan(pilotGrid)).'];
                            % Combine pilot estimates, reference symbol estimates
                            % and those found outwith the resource grid location
                            % and use these to perform 2D interpolation
                            p_est = [p_est outRangeVPs]; %#ok<AGROW>
                        else
                            p_est = [p_est vps]; %#ok<AGROW>
                        end
                    end

                    if (~isempty(p_est))
                        if (strcmpi(cec.InterpType,'None'))
                            % No interpolation, just place pilot estimates 
                            % into the output
                            Htemp = zeros(K,Lsf);
                            p_use_ind = sub2ind(size(Htemp),p_est(1,:),p_est(2,:));
                            Htemp(p_use_ind) = p_est(3,:);
                            Hest(:,1+Lsf*(focalsf-1):Lsf*focalsf,rxANT,tx) = Htemp;
                        else
                            % 2-Dimensional interpolation is carried out to 
                            % estimate the channel between the least 
                            % squares estimates
                            Hest(:,1+Lsf*(focalsf-1):Lsf*focalsf,rxANT,tx) = griddata(p_est(2,:),p_est(1,:),p_est(3,:),1+Lsf*(focalsf-1):Lsf*focalsf,(1:K)',cec.InterpType); %#ok<GRIDD>
                        end
                        
                        % The noise level present can be determined using the noisy
                        % least squares estimates of the channel at pilot symbol
                        % locations and the noise averaged pilot symbol estimates of the
                        % channel.
                        noise = (ls_estimates(3,:)-p_est(3,1:size(ls_estimates,2)));
                        noise = sqrt(scalingVec./(scalingVec+1)).*noise;

                        % Additional averaging in the locations of the 
                        % PUSCH DRS for noise estimation in LTE-A case, to
                        % suppress interference from orthogonal sequences on
                        % other antennas in same time-frequency locations.
                        if (~mod(cec.FreqWindow,12) && (cec.TimeWindow<2))
                            if (size(pusch.PRBSet,1)>1 && strcmpi(cec.Reference,'None')~=1)
                                drsNoise = cell(2,1);
                                for v=1:2
                                    thisNoise=conv(noise(ls_estimates(2,:)==pilotSym(v)),ones(1,12)/sqrt(2));                                                                      
                                    thisNoise(1:12)=[];
                                    thisNoise(end-11:end)=[];
                                    drsNoise{v}=thisNoise;
                                end   
                                noise=(drsNoise{1}-drsNoise{2})/sqrt(2);
                            else
                                % for single PRB case or where PUSCH DRSs 
                                % are not present, this technique is not 
                                % possible, so just rely on the residual 
                                % variance on the averaged pilot estimates.
                                noise=diff(p_est(3,1:size(ls_estimates,2))).*scalingVec(2:end);
                            end
                        end

                        % Taking the variance of the noise present on the pilot symbols
                        % results in a value of the noise power for each transmit and
                        % receive antenna pair
                        noiseVec(rxANT,tx) = mean(noise.*conj(noise));

                    end

                end
            end
        end

        % For channel estimation on the antennas, adjust to give 
        % noise estimates appropriate to the transmitted layers. 
        if (strcmpi(cec.Reference,'Antennas')==1)
            noiseVec=noiseVec/ue.NTxAnts;
        end

        NoisePowerEst = mean(mean(noiseVec));
    else
        NoisePowerEst = NaN;
    end
    
end

% GetPilotEstimates Obtain the least squares estimates of the reference signals
%  [ls_estimates,focalsf,refIndices,pilotSym] = GetPilotEstimates(ue,pusch,tx,rxGrid,cec,refGrid)
%   Extracts the reference signals and calculates their least squares 
%   estimates. The results are placed in a 3xNp matrix containing the 
%   subcarrier and OFDM symbol location, row and column subscripts, and 
%   value. Np is the number of cell specific reference (pilot) symbols per
%   resource grid

%   Copyright 2009-2014 The MathWorks, Inc.

function [ls_estimates,focalsf,refIndices,pilotSym] = GetPilotEstimates(ue,pusch,tx,rxGrid,cec,refGrid)

if (strcmpi(cec.Reference,'Antennas')==1)
    NTx = ue.NTxAnts;
else
    NTx = pusch.NLayers;
end

% Get dimensions of resource grid
Dims = lteULResourceGridSize(ue);
Lsf = Dims(2);

winsizerx = size(rxGrid,2)/Lsf;
if (winsizerx<1)
    error('lte:error','The received grid input must contain at least one subframe.');
end

% Determine window size
winsize = size(refGrid,2)/Lsf;
if (winsize<1)
    error('lte:error','The reference grid input must contain at least one subframe.');
end

if winsize>1
    if isfield(cec,'Window')
        if strcmpi(cec.Window,'Right')
            focalsf = winsize;
            subframeRange = ue.NSubframe-winsize+1:ue.NSubframe;
            if (rem(winsize,1)~=0)
                error('lte:error','The reference grid input must contain a whole number of subframes.');
            end
        elseif strcmpi(cec.Window,'Centred')||strcmpi(cec.Window,'Centered')
            if~mod(winsize,2)
                error('lte:error','Window size must be odd for centred/centered window type.');
            end
            focalsf = floor(winsize/2)+1;
            subframeRange = (ue.NSubframe-floor(winsize/2)):(ue.NSubframe+floor(winsize/2));
        elseif strcmpi(cec.Window,'Left')
            focalsf = 1;
            subframeRange = ue.NSubframe:ue.NSubframe+winsize-1;
        else
            error('lte:error','Window must be Left,Right,Centred/Centered for window sizes>1; see help for details.');
        end
    else
        size3d=@(x)(arrayfun(@(y)size(x,y),1:3));        
        if (all(size3d(rxGrid)==size3d(refGrid)) && all(isnan(refGrid(:))))
            gridstr = 'received';
        else
            gridstr = 'reference';
        end
        error('lte:error',['If the ' gridstr ' grid contains more than one subframe then a channel estimation configuration structure is required including the Window field.']);
    end
else
    focalsf = 1;
    subframeRange = ue.NSubframe;
end

% If locally-generated reference signals are used:
if (strcmpi(cec.Reference,'None')~=1)
    offset = 0;
    drsIndices = [];
    drsSymbols = [];
    for subframe = subframeRange
        ue.NSubframe = mod(subframe,10);
        % Generate linear indices of reference signals
        ind=ltePUSCHDRSIndices(ue,pusch);
        ind=ind(:,tx+1)-(Dims(1)*Lsf*tx);
        drsIndices=[drsIndices ind+(Dims(1)*Lsf*NTx*offset)]; %#ok<AGROW>
        offset = offset +1;
        % Generate expected pilot symbol values
        if (strcmpi(cec.Reference,'Antennas')==1)
            drsTx=ltePUSCHDRS(ue,pusch);
        else
            [~,~,drsTx]=ltePUSCHDRS(ue,pusch);
        end
        drsTx=drsTx(:,tx+1);
        drsSymbols=[drsSymbols drsTx]; %#ok<AGROW>
    end
    [~, row] = ind2sub(size(refGrid),drsIndices);
    pilotSym = [min(row(:,focalsf)) max(row(:,focalsf))];

    % Signal that reference grid in DRS indices are active, even when 
    % DRS symbols might contain zeros.
    refGrid(drsIndices(:)) = 1;
    
else
    
    pilotSym = [];
    
end


% Determine indices of useful symbols in reference grid 
% (non-NaN) and with a reasonable minimum power. The power threshold here
% avoids picking up very noisy REs that may result due to SC-FDMA PUSCH
% precoding - some transmitted REs have inherently low amplitude. 
refGrid(abs(refGrid)<=0.1)=NaN;
refIndices = find(~isnan(refGrid));

% If locally-generated reference signals are used:
if (strcmpi(cec.Reference,'None')~=1)
    % Insert DRS Symbols into correct subframe of refGrid, this may insert
    % zeros for some antennas.
    refGrid(drsIndices(:)) = drsSymbols(:);
end

% Calculate least squares estimates of received data symbols with those
% found in reference grid. Includes handling of producing zeros in 
% locations with zero-valued DRS symbols. 
p_est = rxGrid(refIndices);
ref=refGrid(refIndices);
p_est(ref==0)=0;
ref(ref==0)=1;
p_est=p_est./ref;

% Extract the row and column subscripts of the pilot symbols for entire
% grid
[p_estSC, p_estSym] = ind2sub(size(refGrid),refIndices);

% Reference ls_estimates
ls_estimates = [double(p_estSC(:).') ; double(p_estSym(:).') ; p_est(:).'];

end


%PilotAverage Average reference signals
%   [P_EST] = PilotAverage(CEC,H_EST,LS_EST) performs a moving average of pilot
%   symbols
%
%   LS_EST is a 3xNp matrix containing the least square estimates of the 
%   pilots symbols and their column and row indices within the received
%   grid. 
%   LS_EST = [k;l;p_est]
%   H_EST is an MxN matrix and defines the size of the grid that the
%   averaging will be performed on
%
%   CEC is a structure which defines the type of channel estimation
%   performed. CEC must contain a set of the following fields:
%   FreqWindow      -   Size of window used to average in frequency in
%                       resource elements (odd number, 
%                       or a multiple of 12).
%   TimeWindow      -   Size of window used to average in time in resource
%                       elements (odd number). 
% 
%   The dimensions of the averaging window are defined in structure CEC.
%   The window is defined in terms of Resource Elements, and depending on
%   the size of the averaging window, averaging will be performed in either
%   the time or frequency direction only, or a combination of both creating
%   a square/rectangular window. The pilot to be averaged will always be
%   placed at the center of the window, therefore both FreqWindow and
%   TimeWindow must be odd.

%   Copyright 2009-2014 The MathWorks, Inc.

function [P_EST, scalingVec] = PilotAverage(CEC,H_EST,P_EST)

    if ~mod(CEC.FreqWindow,12) && (CEC.TimeWindow<2)
        
        % Special case where averaging is always applied across (12*X)
        % subcarriers, even at the upper and lower band edges; therefore
        % the first (6*X) symbols at the upper and lower band edge will
        % have the same channel estimate. This operation ensures that
        % averaging is always done on 12 (or a multiple of 12) symbols.
        % This provides the appropriate "despreading" operation required
        % for the case of Uplink MIMO.
        N = CEC.FreqWindow;

        if (~isempty(P_EST))
            
            %Average only SCFDMA symbols which contain DRS symbols. Extract DRS
            %Symbols on a per slot basis and use window of 12 to average.        
            for symbInSlot = unique(P_EST(2,:))

                %Define temporary vector to store the averaged values
                avgVec = zeros(size((P_EST(3,P_EST(2,:)==symbInSlot)).'));
                %Store DRS symbols from relevant slot in vector for easy access
                symbVec = (P_EST(3,P_EST(2,:)==symbInSlot)).';
                %Check length of input at least equal to size of window
                if (length(symbVec) < N)
                    % sprintf ('Input signal must have at least %d elements',N);
                    error('lte:error','Input signal must have at least %d elements.',N);
                end

                %Determines position of window w.r.t element being averaged, and
                %performs the averaging
                for symbNo = 1:length(symbVec)
                    if (symbNo-(N/2+1)<0)
                        avgVec(symbNo)= sum(symbVec(1:N))/N;%numel(symbVec(1:N));
                    elseif (numel(symbVec)<=symbNo+(N/2-1))
                        avgVec(symbNo) = sum(symbVec(end-(N-1):end))/N;%numel(symbVec(end-N:end));
                    else
                        avgVec(symbNo) = sum(symbVec(symbNo-(N/2):symbNo+(N/2-1)))/N;%numel(symbVec(symbNo-(N/2):symbNo+(N/2-1)));
                    end
                end

                %Update P_EST with averaged values
                P_EST(3,P_EST(2,:)==symbInSlot)= avgVec;
            end
            
        end
        
        %This vector is used to scale the noise by the number of averaging
        %elements in the window
        scalingVec = ones(size(P_EST(3,:)))*N;        
        
    else
        
        % Check window dimensions are odd placing the averaged pilot at the
        % center
        if ~mod(CEC.FreqWindow,2)||~mod(CEC.TimeWindow,2)
            error('lte:error','For rectangular pilot averaging, both FreqWindow and TimeWindow must be odd.');
        end
        
        % Define an empty resource grid
        reGrid = zeros(size(H_EST,1),size(H_EST,2));
        % Place the pilot symbols back into the received grid
        reGrid(sub2ind(size(reGrid),P_EST(1,:),P_EST(2,:))) = P_EST(3,:);
        % Define convolution window
        kernel = ones(CEC.FreqWindow,CEC.TimeWindow);
        % Average by performing 2-dimensional convolution
        reGrid = conv2(reGrid,kernel,'same');
        % Extract only pilot symbol location values and set the rest of
        % the grid to zero
        tempGrid = zeros(size(reGrid));
        tempGrid(sub2ind(size(reGrid),P_EST(1,:),P_EST(2,:))) = reGrid(sub2ind(size(reGrid),P_EST(1,:),P_EST(2,:)));
        reGrid = tempGrid;
        % Normalize pilot symbol values after convolution
        [reGrid, scalingVec] = normalisePilotAverage(CEC,P_EST,reGrid);
        % Place averaged values back into pilot symbol matrix
        P_EST(3,:) = reGrid(sub2ind(size(tempGrid),P_EST(1,:),P_EST(2,:)));
        
    end
    
end

function [avgGrid, scalingVec] = normalisePilotAverage(CEC,p_est,reGrid)

% Determine total number of pilots within half a subframe
nPilots = length(p_est);

avgGrid = zeros(size(reGrid));
scalingVec = zeros(1,size(p_est,2));

for n = 1:nPilots
    % Determine in which subcarrier and OFDM symbol pilot is located 
    sc = p_est(1,n);
    sym = p_est(2,n);
    % Determine number of REs to look at either side of pilot 
    % symbol in time and frequency
    half_freq_window = floor(CEC.FreqWindow/2);
    half_time_window = floor(CEC.TimeWindow/2);
    % Set the location of the window at the back of the pilot 
    % to be averaged in frequency direction
    upperSC = sc-half_freq_window;
    % If this location is outwith the grid dimensions set it to
    % lowest subcarrier value
    if upperSC<1
        upperSC = 1;
    end
    % Set the location of the window in front of the pilot to 
    % be averaged in frequency direction
    lowerSC =  sc+half_freq_window;
    % If this location is outwith the grid dimensions set it to
    % highest subcarrier value
    if lowerSC>size(reGrid,1)
        lowerSC = size(reGrid,1);
    end
    % Set the location of the window in front of the pilot to 
    % be averaged in time direction
    leftSYM = sym-half_time_window;
    % If this location is outwith the grid dimensions set it to
    % lowest OFDM symbol value
    if leftSYM<1
        leftSYM = 1;
    end
    % Set the location of the window at the back of the pilot 
    % to be averaged in time direction
    rightSYM = sym+half_time_window;
    % If this location is outwith the grid dimensions set it to
    % highest OFDM symbol value
    if rightSYM>size(reGrid,2)
        rightSYM = size(reGrid,2);
    end

    % Define the window to average using the determined 
    % subcarrier and OFDM symbol values
    avgVec = reGrid(upperSC:lowerSC,leftSYM:rightSYM);

    % Remove the zero values from the window so that the average is
    % calculated using only valid pilots
    avgVec = avgVec(avgVec~=0);
    % Average the desired pilot using all pilots within 
    % averaging window
    if isempty(avgVec)
        avgVec = 0;
    end

    avgGrid(sc,sym) = reGrid(sc,sym)/numel(avgVec);
    scalingVec(n) = numel(avgVec);
end
end

%CreateVirtualPilots 
% Creates virtual pilots on the edges of the resource grid to improve
% interpolation results
function vps = CreateVirtualPilot(ue,p_est,focalsf)

    % Determine how many symbols are in current subframe and add 2 to
    % account for edges
    dims=lteULResourceGridSize(ue);
    Lsf = dims(2);
    symbols=dims(2)/2;
    totSym = symbols*2+2;
    % Total number of subcarriers with reference signals
    maxSC = max(p_est(1,:));
    % Determine total number of pilot symbol carrying subcarriers
    totSC = (length(p_est(2,:)==p_est(2,1)))/2;
    scOffset = min(p_est(1,:));
     
    pilotsym = unique(p_est(2,:));
    vps = [];
    for i = 1:totSC
        newPil = interp1(p_est(2,i)+1:symbols:p_est(2,i+totSC)+1,[p_est(3,i) p_est(3,i+totSC)],1+(Lsf*(focalsf-1)):(Lsf*(focalsf-1))+totSym,'linear','extrap');
        vps = [vps [i+scOffset-1;(Lsf*(focalsf-1));newPil(1)] [i+scOffset-1;Lsf*(focalsf-1)+totSym-1;newPil(end)]]; %#ok<AGROW>
    end
    for j = pilotsym
        pil = interp1(min(p_est(1,:))+1:max(p_est(1,:))+1,p_est(3,p_est(2,:)==j) ,min(p_est(1,:)):max(p_est(1,:))+2,'nearest','extrap');
        vps = [vps [scOffset-1;j;pil(1)] [maxSC+1;j;pil(end)]]; %#ok<AGROW>
    end
end

%CreateEdgeVirtualPilots 
% Creates virtual pilots on the edges of the resource grid to improve
% interpolation results. Also adds virtual pilots around the current pilot
% estimates.
function vps = createEdgeVirtualPilots(ue,p_est)
       
    % Determine dimensions of current subframe
    dims=lteULResourceGridSize(ue);
    K=dims(1);
    L=dims(2);

    % Calculate virtual pilots beyond upper and lower bandwidth edge based
    % on pilots present on subcarriers. Also create virtual pilots 1/2 RB
    % above and below extent of current pilots. 
    vps=createDimensionVirtualPilots(p_est,2,unique([-6 min(p_est(1,:)-6) max(p_est(1,:)+6) K+6]));                    

    % Combine these frequency-direction VPs with original pilots.
    temp = [p_est vps];

    % Calculate virtual pilots beyond start and end of subframe based on
    % pilots and frequency-direction VPs present in symbols. Also create 
    % virtual 1 OFDM symbol before and after extent of current pilots. 
    vps = [vps createDimensionVirtualPilots(temp,1,unique([-1 min(p_est(2,:)-1) max(p_est(2,:)+1) L+1]))];
    
end

% dim=1 adds VPs in time
% dim=2 adds VPs in frequency
function vps = createDimensionVirtualPilots(p_est,dim,points)

    vps=[];        
    
    pilots = unique(p_est(dim,:));        
    for i = pilots
        x=find(p_est(dim,:)==i);
        rep=1+double(length(x)==1);
        adj=(1:rep)-rep;
        pil = interp1(p_est(3-dim,x)+adj,repmat(p_est(3,x),1,rep),points,'linear','extrap');
        for j=1:length(points)
            vps = [ vps [i;points(j);pil(j)] ];         %#ok<AGROW>
        end
    end
    
    if (dim==2)
        vps = vps([2 1 3],:);
    end
    
end

function defaultValueWarning(field,value)
    s=warning('query','backtrace');
    warning off backtrace;        
    warning('lte:defaultValue','Using default value for parameter field %s (%s)',field,value);
    warning(s); 
end
