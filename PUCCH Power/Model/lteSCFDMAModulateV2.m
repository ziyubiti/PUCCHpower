
function [timeDomainSig,infoScfdma] = lteSCFDMAModulateV2(rmc,frame)

[timeDomainSig,infoScfdma] = lteSCFDMAModulate(rmc,frame);
timeDomainSig = timeDomainSig*sqrt(infoScfdma.Nfft);