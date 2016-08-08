%lteSCFDMADemodulate SC-FDMA demodulation
%   GRID = lteSCFDMADemodulate(UE,WAVEFORM) performs SC-FDMA demodulation
%   of the time domain waveform WAVEFORM given ue-specific settings
%   structure UE which must include the following fields:
%   NULRB          - Number of uplink resource blocks
%   CyclicPrefixUL - Optional. Cyclic prefix length 
%                    ('Normal'(default),'Extended')
%
%   The demodulation performs one FFT operation per received SC-FDMA
%   symbol, to recover the received subcarrier values which are then used
%   to construct each column of the output resource array GRID. The FFT is
%   positioned part-way through the cyclic prefix, to allow for a certain
%   degree of channel delay spread whilst avoiding the overlap between
%   adjacent OFDM symbols. The input the FFT is also shifted by half of one
%   subcarrier. The particular position of the FFT chosen here avoids the
%   SC-FDMA symbol overlapping used in the <a href="matlab:
%   help('lteSCFDMAModulate')">lteSCFDMAModulate</a> function. 
%   Given that the FFT is performed away from the original zero phase point
%   on the transmitted subcarriers, a phase correction is applied to each
%   subcarrier after the FFT.
%
%   Note that the sampling rate of the time domain waveform WAVEFORM must
%   be the same as used in the lteSCFDMAModulate modulator function for the
%   specified number of resource blocks NULRB. WAVEFORM must also be
%   time-aligned such that the first sample is the first sample of the
%   cyclic prefix of the first SC-FDMA symbol in a subframe.
%
%   GRID = lteSCFDMADemodulate(UE,WAVEFORM,CPFRACTION) allows the
%   specification of the position of the demodulation through the cyclic
%   prefix. CPFRACTION is between 0 and 1, with 0 representing the start of
%   the cyclic prefix and 1 representing the end of the cyclic prefix. The
%   default value is 0.55, which allows for the default level of windowing
%   in the <a href="matlab:
%   help('lteSCFDMAModulate')">lteSCFDMAModulate</a> function.
%
%   Example: 
%   SC-FDMA demodulation of Uplink FRC A3-2.
%   
%   frc = lteRMCUL('A3-2');
%   waveform = lteRMCULTool(frc,randi([0,1],frc.PUSCH.TrBlkSizes(1),1));   
%   reGrid = lteSCFDMADemodulate(frc,waveform);
%
%   See also lteSCFDMAModulate, lteSCFDMAInfo, lteULFrameOffset,
%   lteULFrameOffsetPUCCH1, lteULFrameOffsetPUCCH2, lteULFrameOffsetPUCCH3,
%   lteULChannelEstimate, lteULPerfectChannelEstimate,
%   lteULChannelEstimatePUCCH1, lteULChannelEstimatePUCCH2,
%   lteULChannelEstimatePUCCH3.

%   Copyright 2010-2014 The MathWorks, Inc.

function reGrid = lteSCFDMADemodulate(ue,waveform,varargin)

    if (nargin==3)
      cpFraction = varargin{1};
    else
      cpFraction = 0.55;
    end
    
    % Compute the number of samples per subframe and FFT size from
    % the UE settings.
    info=lteSCFDMAInfo(ue);
    samplesPerSubframe=info.SamplingRate*0.001;
    nFFT=double(info.Nfft);
    
    % Get the cyclic prefix lengths for one slot for the cyclic prefix
    % type specified in UE.
    cpLengths=double(info.CyclicPrefixLengths);    
    cpLengths=cpLengths(1:length(cpLengths)/2);
    
    % Calculate the number of symbols in a slot from the cyclic prefix
    % lengths, and the number of whole subframes in WAVEFORM.
    nSymbols=length(cpLengths);    
    nSubframes=fix(size(waveform,1)/samplesPerSubframe);        

    % Calculate position of the first active subcarrier in the FFT output,
    % and the total number of active subcarriers, according to the FFT
    % size and number of resource blocks. 
    firstActiveSC=(nFFT/2)-(ue.NULRB*6)+1;
    totalActiveSC=ue.NULRB*12;
    
    % Create an empty output GRID of the correct size.
    dims=lteULResourceGridSize(ue,size(waveform,2));
    dims(2)=dims(2)*nSubframes;    
    reGrid=zeros(dims);    
    idx=0:nFFT-1;   
    
    % For each subframe in WAVEFORM, each slot in a subframe and each
    % symbol in a slot:  
    i=1;    
    offset=0;
    for subframe=0:nSubframes-1 
        for slot=0:1
            for symbol=0:nSymbols-1                
                                
                % Get cyclic prefix length in samples for the current symbol.
                cpLength=cpLengths(symbol+1);
                
                % Position the FFT part of the way through the cyclic
                % prefix; the value of cpFraction should ensure that the 
                % reception is unaffected by the windowing applied in the 
                % lteOFDMModulate function. Default is 55%.                
                fftStart=fix(cpLength*cpFraction);
                
                % Create vector of phase corrections, one per FFT sample,
                % to compensate for FFT performed away from zero phase
                % point on the original subcarriers.
                phaseCorrection=repmat(exp(-1i*2*pi*(cpLength-fftStart)/nFFT*idx)',1,size(waveform,2));               
                
                % Extract the appropriate section of WAVEFORM, perform the
                % FFT and half-subcarrier shifting and apply the phase correction.
                halfsc=repmat(exp(1i*pi/nFFT*(idx+fftStart-cpLength))',1,size(waveform,2));                
                fftOutput=fftshift(fft(waveform(offset+fftStart+(1:nFFT),:).*halfsc).*phaseCorrection,1)./sqrt(nFFT);
                
                % Extract the active subcarriers for each antenna from the 
                % FFT output.
                activeSCs=fftOutput(firstActiveSC:firstActiveSC+totalActiveSC-1,:);
                
                % Assign the active subcarriers into the appropriate column
                % of the received GRID, for each antenna.
                reGrid(SCFDMASymbolIndices(reGrid,i))=activeSCs;
                
                % update counter of overall symbol number and position in
                % the input WAVEFORM.
                i=i+1;
                offset=offset+nFFT+cpLength;
            end
        end
    end    
end

%SCFDMASymbolIndices Generates indices for a given SCFDMA symbol of a resource matrix.
%   IND = SCFDMASymbolIndices(GRID,SYMBOL) gives the indices for SC-FDMA symbol number
%   SYMBOL of the resource array GRID, in per-antenna linear indices form. 
function ind = SCFDMASymbolIndices(reGrid,symbol)
    nSCs = size(reGrid,1);
    nAnts = size(reGrid,3);
    firstSCs = sub2ind(size(reGrid),ones(1,nAnts),repmat(symbol,1,nAnts),1:nAnts);
    ind = repmat(firstSCs,nSCs,1)+repmat((0:nSCs-1).',1,nAnts);    
end

