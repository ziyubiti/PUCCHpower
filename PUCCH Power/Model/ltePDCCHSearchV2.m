%ltePDCCHSearch PDCCH downlink control information search
%   [DCISTR,DCIBITS] = ltePDCCHSearch(ENB,CHS,SOFTBITS) recovers DCI
%   message structures DCISTR and corresponding vectors of DCI message bits
%   DCIBITS after blind decoding the multiplexed PDCCHs within the control
%   region given by input vector of soft bits SOFTBITS, cell-wide settings
%   ENB and User Equipment (UE) related configuration structure CHS.
%   
%   ENB must be a structure including the fields:
%   NDLRB      - Number of downlink resource blocks (see note below)
%   NULRB      - Number of uplink resource blocks (see note below)
%   NSubframe  - Subframe number
%   CellRefP   - Number of cell-specific reference signal antenna ports 
%                (1,2,4)
%   DuplexMode - Optional. Frame structure ('FDD'(default),'TDD')
%   
%   CHS must be a structure including the field:
%	RNTI - RNTI value (16 bits)
%   
%   Note that at least one of the NDLRB and NULRB parameter pair must be
%   provided. If either is not provided then the function assumes symmetric
%   link bandwidth where NDLRB and NULRB are equal. No warning message will
%   be displayed in this event.
%   
%   The UE is required to monitor multiple PDCCHs within the control
%   region. The UE is only informed of the width (in OFDM symbols) of the
%   control region within a subframe and is not aware of the exact location
%   of PDCCHs relevant to it. The UE finds the PDCCHs relevant to it by
%   monitoring a set of PDCCH candidates (a set of consecutive control
%   candidate elements (CCEs) on which PDCCH could be mapped) in every
%   subframe (see <a href="matlab: 
%   help('ltePDCCHSpace')">ltePDCCHSpace</a> for details). This is referred
%   to as blind decoding.
%   
%   To simplify the decoding task at the UE, the whole control region is
%   sub-divided into common and UE-specific search spaces which the UE
%   monitors (monitor implies attempting to decode each PDCCH). Each search
%   space comprises 2, 4 or 6 PDCCH candidates whose data length depends on
%   its corresponding PDCCH format; each PDCCH must be transmitted on 1, 2,
%   4 or 8 CCE (1 CCE = 72 bits). The common search space is limited to
%   only two aggregation levels, 4 and 8, while the UE-specific search
%   space can have an aggregation level of 1, 2, 4, or 8.
%   
%   The common search space carries control information common to all UEs
%   and is therefore monitored by all UEs within a cell. The common control
%   information carries initial important information including paging
%   information, system information and random access procedures. The UE
%   monitors the common search space by demasking each PDCCH candidate
%   with different RNTIs e.g. P-RNTI, SI-RNTI, RA-RNTI etc.
%   
%   In the UE-specific search space the UE finds the PDCCH relevant to it
%   by monitoring a set of PDCCH candidates in every subframe. If no CRC
%   error is detected when the UE demasks a PDCCH candidate with its RNTI
%   (16-bit C-RNTI value), the UE determines that the PDCCH candidate
%   carries its own control information.
%   
%   The number and location of candidates within a search space is
%   different for each PDCCH format. There are four PDCCH formats i.e. 0,
%   1, 2 or 3. If the UE fails to decode any PDCCH candidates for a given
%   PDCCH format then it tries to decode candidates for another PDCCH
%   format.
%   
%   DCI messages are blindly decoded on the basis of their lengths. The
%   lengths and order in which they are searched for is provided by 
%   <a href="matlab:help('lteDCIInfo')">lteDCIInfo</a>. If one or more messages have the same length then the
%   message format that is first in the list is used to decode the message.
%   The other potential message formats are ignored. The transmission mode
%   (TM) is not taken into account during blind search and no DCI message
%   format is filtered on the basis of transmission mode. Format 3 and 3A
%   (power adjustment commands for PUSCH and PUCCH) are not searched by
%   this function.
%   
%   DCISTR is a cell array of structures containing the fields associated
%   with one or more decoded DCI message(s). As multiple PDCCHs can be
%   transmitted in a subframe the UE has to monitor all possible PDCCHs
%   directed at it. If more than one PDCCH is directed to the UE or
%   successfully decoded then DCISTR will contain that number of decoded
%   DCI messages.
%   
%   DCIBITS is a cell array containing one or more vectors of bit values
%   corresponding to successfully decoded DCI messages.
%   
%   Example:
%   The PDCCH symbols are extracted from the transmit resource grid and
%   decoded. The UE monitors the common search space by demasking the PDCCH
%   candidate with the configured RNTI. The DCI cell array rxDCI and
%   corresponding cell array of DCI message bit vectors rxDCIBits are
%   recovered after blind decoding the multiplexed PDCCH within the control
%   region.
%   
%   rmc = lteRMCDL('R.0');            
%   [~,txGrid] = lteRMCDLTool(rmc,[1;0;0;1]);
%   % extract and decode PDCCH bits
%   pdcchSymbols = txGrid(ltePDCCHIndices(rmc));
%   rxPdcchBits = ltePDCCHDecode(rmc,pdcchSymbols); 
%   % PDCCH blind search, demask PDCCH candidate using RNTI
%   ueConfig.RNTI = rmc.PDSCH.RNTI; 
%   [rxDCI,rxDCIBits] = ltePDCCHSearch(rmc,ueConfig,rxPdcchBits);
%   % extract and display first DCI message structure
%   decDCI = rxDCI{1}
%   
%   See also ltePDCCH, ltePDCCHDecode, ltePDCCHIndices, ltePDCCHInterleave,
%   ltePDCCHDeinterleave, ltePDCCHInfo, ltePDCCHSpace, ltePDCCHPRBS.

%   Copyright 2009-2014 The MathWorks, Inc.

function [decDCI,decDCIBits,dciCfg] = ltePDCCHSearchV2(enbConfig,ueConfig,inBits)
        
    % DCI configuration structure
    dciConfig = enbConfig;
    
    % Checking mandatory parameter
    if(~isfield(enbConfig,'CellRefP'))
        error('lte:error','The function call (ltePDCCHSearch) resulted in an error: Could not find a structure field called CellRefP.');
    end
    if(~isfield(ueConfig,'RNTI'))
        error('lte:error','The function call (ltePDCCHSearch) resulted in an error: Could not find a structure field called RNTI.');
    end
    
    % Input size validation
    if(~(length(inBits)>=72))
        error('lte:error','Input soft bits should at least be a CCE in length (72 bits).');
    end
    
    % Deducing total number of REGs (1 REG = 4 REs = 8 bits)
    enbConfig.NREG = floor(length(inBits)/(72/9));
    
    % DCI formats for Common and UE-Specific search space
    % Common search space DCI formats
    dciFormats{1}={'Format0'; 'Format1A'; 'Format1C'};
    
    % UE-Specific search space DCI formats
    dciinfo = lteDCIInfo(enbConfig);
    dcimessages = fieldnames(dciinfo);
    uespecific_exclusionList = {'Format3' 'Format3A'};
    dciFormats{2} = setdiff(dcimessages,uespecific_exclusionList);
    
    % PDCCH format for Common search space can either be 2 or 3 (i.e.
    % aggregation level of 4 or 8 CCEs
    startingPdcchFormat = 2;
    
    % Intermediate helping variables
    idx=1;
    decDCI={};
    decDCIBits={};
    reservedLoc=[];
    dciCfg = struct('PDCCHFormat', NaN);       % initial default
    
    if (ueConfig.RNTI == 65535) || (ueConfig.RNTI == 65534)
        endType = 1;        
    else
        endType = 2;
    end;
    
    for searchType=1:endType
        % UE-specific search space
        if(searchType == 2)
            pdcchConfig.RNTI = ueConfig.RNTI;
            
            % PDCCH format for ue-specific search space can be 0,1,2 or 3
            startingPdcchFormat = 0;
        end;
        startingPdcchFormat = 2; % unit test
        for pdcchFormat = startingPdcchFormat:1:3
            pdcchConfig.PDCCHFormat = pdcchFormat;

            % Performs common and/or ue-specific search
            pdcchCandidates = ltePDCCHSpace(enbConfig,pdcchConfig,{'bits','1based'});

            % PDCCH candidates need not to be unique so picking the unique set
            % of candidates to optimize the search
            pdcchCandidates = unique(pdcchCandidates,'rows');

            pdcchCandidatesDims = size(pdcchCandidates);
            noOfCandidates = pdcchCandidatesDims(1);

            for candidate=1:noOfCandidates
                if(sum(reservedLoc == pdcchCandidates(candidate,1)/72) == 0)
                    if((pdcchCandidates(candidate,1)<length(inBits)) &&  (pdcchCandidates(candidate,2)<=length(inBits)))
                        input = inBits(pdcchCandidates(candidate,1):pdcchCandidates(candidate,2));
                        for dciFormatIdx=1:length(dciFormats{searchType}) % Iterating through all DCI formats
                            dciConfig.DCIFormat=dciFormats{searchType}{dciFormatIdx};
                            [dciMessageBits,decRnti] = lteDCIDecode(dciConfig,input);
                            if(ueConfig.RNTI == decRnti && ~any(reservedLoc == pdcchCandidates(candidate,1)/72))
                                % Creating DCI message for decoded bits
                                dciCfg.PDCCHFormat = pdcchFormat;  % level 1 2 4 8 for 0 1 2 3
                                dciCfg.DCIFormat = dciConfig.DCIFormat;
                                dciCfg.candidate = candidate;
                                [dciMessage,dciMessageBits] = lteDCI(enbConfig,dciMessageBits);
                                decDCIBits{idx} = dciMessageBits;
                                decDCI{idx} = dciMessage;
                                 reservedLoc = [reservedLoc,(pdcchCandidates(candidate,1)/72:(pdcchCandidates(candidate,1)/72)+(2^pdcchFormat)-1)]; %#ok<AGROW>
                                idx = idx + 1;
                                break;
                            end
                        end
                    end
                end
            end
        end
    end
end
