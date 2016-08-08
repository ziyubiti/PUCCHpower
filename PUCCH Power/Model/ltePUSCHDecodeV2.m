%ltePUSCHDecode Physical uplink shared channel decoding
%   [CWS,SYMBOLS] = ltePUSCHDecode(...) returns a soft bit vector, or cell
%   array of soft bit vectors, CWS containing the received codeword
%   estimates and received constellation of complex symbol vector SYMBOLS
%   resulting from performing the inverse of the Physical Uplink Shared
%   Channel (PUSCH) processing (TS 36.211 5.3, for details see <a
%   href="matlab:help('ltePUSCH')">ltePUSCH</a>).
%   Codeword or codewords CWS are optionally scaled by channel state
%   information (CSI) calculated during the equalization process.
%
%   [CWS,SYMBOLS] = ltePUSCHDecode(UE,CHS,SYM) returns a soft bit vector,
%   or cell array of soft bit vectors, CWS containing the received codeword
%   estimates and received constellation of complex symbol vector SYMBOLS
%   resulting from decoding of the Physical Uplink Shared Channel (PUSCH)
%   complex symbols SYM for UE-specific settings UE and channel
%   transmission configuration structure or structure array CHS.
%   
%   UE must be a structure including the fields:
%   NCellID        - Physical layer cell identity
%   NSubframe      - Subframe number
%   RNTI           - Radio Network Temporary Identifier (16-bit)
%   CyclicPrefixUL - Optional. Cyclic prefix length 
%                    ('Normal'(default),'Extended') 
%   NTxAnts        - Optional. Number of transmission antennas 
%                    (1(default),2,4).
%   Shortened      - Optional. If 1, the last symbol of the subframe is not
%                    used and rate matching is adjusted accordingly This is
%                    required for subframes with possible SRS transmission.
%                    Default 0.
%   
%   CHS is the PUSCH channel-specific structure including the fields:
%   Modulation - Modulation format ('QPSK','16QAM','64QAM')
%   PRBSet     - A 1- or 2-column matrix, containing the Physical Resource 
%                Block (PRB) indices corresponding to the slot wise
%                resource allocations for this PUSCH.
%   NLayers    - Optional. Number of transmission layers 
%                (total or per codeword) (1(default),2,3,4)
%   Only required for NTxAnts=2 or NTxAnts=4:
%      PMI     - Optional. Scalar Precoder Matrix Indication to be used  
%                during precoding. (0(default)...23), depending on NTxAnts 
%                and NLayers. (See <a
%                href="matlab:help('lteULPMIInfo')">lteULPMIInfo</a>)
%
%   If UCI is present in the transmitted PUSCH to be decoded, the following 
%   optional fields should be configured in the CHS structure:
%   ORI   - Optional. Number of uncoded RI bits (default 0)
%   OACK  - Optional. Number of uncoded HARQ-ACK bits (default 0)
%   QdRI  - Optional. Number of coded RI symbols in UL-SCH (Q'_RI) 
%           (default 0)
%   QdACK - Optional. Number of coded HARQ-ACK symbols in UL-SCH  (Q'_ACK) 
%           (default 0)
%
%   Multiple codewords can be parameterized by two different forms of the
%   CHS structure. Each codeword can be defined by separate elements of a
%   (1-by-2) structure array, or the codeword parameters can be combined
%   together in the fields of a single scalar (1-by-1) structure (any
%   scalar field values apply to both codewords and a scalar NLayers is the
%   total number). For further details, see <a href="matlab: 
%   web([docroot '/lte/gs/ul-sch-parameterization.html'])">UL-SCH Parameterization</a>.
%   
%   If UCI control information (RI or HARQ-ACK) is present in the received
%   complex PUSCH symbols, then this function performs the descrambling of
%   the placeholder bits by establishing the correct locations with the
%   help of the UCI related parameters present in CHS.
%
%   For PRBSet, if a column vector is provided, the resource allocation is
%   the same in both slots of the subframe; the 2-column matrix can be used
%   to specify differing PRBs for each slot in a subframe. Note that the
%   PRB indices are 0-based.
%
%   SYM is an M-by-P matrix or M-by-NU matrix where M is the number of
%   symbols per antenna or layer, P is the number of transmission antennas
%   (NTxAnts) and NU is the number of transmission layers (NLayers). For a
%   single-antenna transmission (NTxAnts=1), P=NU=1 and SYM must be M-by-1
%   and contain the single-antenna PUSCH symbols for decoding. For P>1,
%   when SYM is M-by-P, decoding is performed using pseudo-inverse based
%   deprecoding for spatial multiplexing. For P>1, when SYM is M-by-NU,
%   decoding is performed without deprecoding; the input SYM is assumed to
%   already be deprecoded for example by having performed channel
%   estimation against the transmit layer DRS sequences and performing
%   equalization of the received symbols using that channel estimate to
%   yield SYM. The case where P>1 and P=NU is ambiguous as to whether or
%   not deprecoding is required here; this function does apply deprecoding
%   in that case.
%
%   [CWS,SYMBOLS] = ltePUSCHDecode(UE,CHS,SYM,HEST,NOISEEST) performs the
%   decoding of the complex PUSCH symbols SYM using UE-specific settings
%   UE, channel transmission configuration CHS, the channel estimate HEST
%   and the noise estimate NOISEEST. In this case SYM is an M-by-NRxAnts
%   matrix where M is the number of symbols per antenna and NRxAnts is the
%   number of receive antennas. For UE.NTxAnts>1, the reception is
%   performed using an MMSE equalizer, equalizing between transmitted and
%   received layers. For UE.NTxAnts=1, the reception is performed using
%   MMSE equalization on the received antennas.
%      
%   HEST is a 3-dimensional M-by-NRxAnts-by-NTxAnts array where M is the
%   number of symbols per antenna, NRxAnts is the number of receive
%   antennas, and NTxAnts is the number of transmit antennas ports, given
%   by UE.NTxAnts.
%   
%   NOISEEST is an estimate of the noise power spectral density per RE on
%   received subframe; such an estimate is provided by the function
%   <a href="matlab:help('lteULChannelEstimate')">lteULChannelEstimate</a>.
%
%   [CWS,SYMBOLS] = ltePUSCHDecode(UE,CHS,SYM,HEST,NOISEEST,ALG) is the
%   same as above except it provides control over weighting the output soft
%   bits BITS with Channel State Information (CSI) calculated during the
%   equalization stage using algorithmic configuration structure ALG.
%   
%   ALG must be a structure including the field:
%   CSI - Optional. Determines if soft bits should be weighted by CSI 
%         ('Off','On'(default)).
%   
%   Example:
%   Decode the PUSCH modulation symbols contained in the output of a
%   Fixed Reference Channel (FRC).
%   
%   frc = lteRMCUL('A3-2');
%   trData = randi([0,1],frc.PUSCH.TrBlkSizes(1),1);
%   [waveform,reGrid] = lteRMCULTool(frc,trData);
%   puschIndices = ltePUSCHIndices(frc,frc.PUSCH);
%   rxCw = ltePUSCHDecode(frc,frc.PUSCH,reGrid(puschIndices));
%
%   See also ltePUSCH, ltePUSCHIndices, ltePUSCHDRS, ltePUSCHDRSIndices,
%   ltePUSCHDeprecode, lteULDescramble, lteULDeprecode, lteULSCHDecode.

%   Copyright 2010-2014 The MathWorks, Inc.

function [cws,delayered] = ltePUSCHDecodeV2(ue,chsin,sym,varargin)
    
    % expand 'chsin' into a structure array if necessary.    
    opfields = {'NLayers','ORI','OACK','QdRI','QdACK'};
    chsin = saLteParamExpand(chsin,['Modulation',opfields],opfields,{1,0,0,0,0,0});
        
    % if 'chsin' is (now) a structure array:
    if (numel(chsin)>1)
        % manipulate 'chsin' structure array to create single 'chs'
        % compatible with all steps except lteULDescramble.
        chs=chsin(1);
        if (isfield(chs,'NLayers'))
            chs.NLayers=sum([chsin.NLayers]);            
        end            
        chs.Modulation={chsin.Modulation};
    else               
        % use 'chsin' directly as 'chs'. 
        chs=chsin;        
    end

    % Provide default value for CSI field if absent
    if(nargin>5)
        alg = varargin{3};
        if (~isfield(alg,'CSI'))
            alg.CSI='On';
        end
    else
        alg.CSI='On';
    end
    
    if (~isfield(ue,'NTxAnts'))
        ue.NTxAnts=1;
        warning('lte:defaultValue','Using default value for parameter field NTxAnts (1)');
    end  
    
    if(nargin == 3)
        
        % Perform spatial deprecoding of received symbols, if appropriate.
        P=size(sym,2);
        if (P==ue.NTxAnts)
            layers = ltePUSCHDeprecode(chs,sym);       
        else
            layers = sym;
        end
               
    elseif (nargin >= 5)
       
        hest=varargin{1};
        noiseEst=varargin{2}; 
        
        if (ue.NTxAnts>1)
            % MIMO equalization and deprecoding (MMSE based).
            [layers,csi] = lteEqualizeULMIMO(ue,chs,sym,hest,noiseEst);
        else           
            % MMSE for single antenna.
            [puschRx,csi] = lteEqualizeMMSE(sym,hest,noiseEst);               
            layers = ltePUSCHDeprecode(chs,puschRx);
        end
        
    end
    
    
    % Perform SC-FDMA deprecoding of layers 
    NPRB = size(chs.PRBSet,1);
    deprecoded = lteULDeprecode(layers,NPRB);
    if (nargin == 3)
        csi = ones(size(deprecoded));
    else
        % Form CSI by deprecoding the CSI, normalising for the effect of
        % the integration in the IDFT, taking the DC component (every Kth
        % value, representing the average CSI over those K values) for each
        % OFDM symbol and layer, repeating these values for all K
        % subcarriers in each OFDM symbol and each layer, and re-arranging
        % back into multiple columns for the layers.
        K = NPRB*saLteNscRB();
        csi = lteULDeprecode(csi,NPRB)/sqrt(K);        
        csi = reshape(repmat(csi(1:K:end).',1,K).',size(deprecoded));
    end
    
    % Perform layer demapping of the SC-FDMA deprecoded symbols 
    delayered = lteLayerDemap(chs,deprecoded);
    csi = lteLayerDemap(chs,csi);
        
    % Perform demodulation of the deprecoded symbols    
    if(iscell(chs.Modulation))
        ncw = length(chs.Modulation);
        demod = cell(1,ncw);
        for cwIdx=1:ncw            
            demod{cwIdx} = lteSymbolDemodulate(delayered{cwIdx},chs.Modulation{cwIdx});
            if (~isempty(demod{cwIdx}))
                Qm = length(demod{cwIdx})/length(delayered{cwIdx});
                csi{cwIdx} = repmat(csi{cwIdx}.',Qm,1);
                csi{cwIdx} = reshape(csi{cwIdx},numel(csi{cwIdx}),1);
            end
        end
    else        
        delayered=delayered{1};
        csi=csi{1};
        % bj:
        if 0
            figure();
            scatter(real(delayered),imag(delayered));
            grid on;
        end;
        
        demod = lteSymbolDemodulate(delayered,chs.Modulation);
        if (~isempty(demod))
            Qm = length(demod)/length(delayered);
            csi = repmat(csi.',Qm,1);
            csi = reshape(csi,numel(csi),1);        
        end
    end
    
    % descramble to produce the codeword estimate
    cws = lteULDescramble(ue,chsin,demod);
    
    % Scaling LLRs by CSI
    if (strcmpi(alg.CSI,'On')==1)
        if (iscell(cws))
            ncw = numel(cws);
            for cwIdx=1:ncw
                cws{cwIdx}=cws{cwIdx}.*csi{cwIdx};
            end
        else
            cws=cws.*csi;
        end
    end
    
end
