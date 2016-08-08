% coarse freq offset estimation,[-7.5K  7.5K]
% flag = 1, UL, the cp part and the last symbole are phase inversed ,
% because the 7.5K half sc freq offset

function [freqoffset_cp] = FreOffsetEstimateCp(timeDomainSig,infoScfdma,ULflag)

nFFT = infoScfdma.Nfft;
CP_len = infoScfdma.CyclicPrefixLengths;

q = size(timeDomainSig,2);


for i = 0:length(CP_len)/2-1
    calc_len(i+1) = CP_len(1,i+1)/2;
    if i==0      
        start_index(i+1) = calc_len(i+1);
    else
        start_index(i+1) = calc_len(i+1) + CP_len(1,1) + (i-1)*CP_len(1,2)+ nFFT*i;
    end;
end;

TS1_start_index = start_index + size(timeDomainSig,1)/2;
start_index_sf = [start_index TS1_start_index];
calc_len_sf = [calc_len calc_len];


data1 = [];
data2 = [];
% data3 = [];
% data4 = [];
for j = 1:q
    for k = 1:length(CP_len)
        tmp1 = timeDomainSig(start_index_sf(1,k)+1:start_index_sf(1,k)+calc_len_sf(k),j);
        tmp2 = timeDomainSig(start_index_sf(1,k)+1+nFFT:start_index_sf(1,k)+calc_len_sf(k)+nFFT,j);
        data1 = [data1;tmp1];
        data2 = [data2;tmp2];
    end;
%     data3 = [data3 data1];
%     data4 = [data4 data2]; 
end;

if ULflag
    xcorr = conj(-data1).*data2;
else
    xcorr = conj(data1).*data2;
end;

phase_xcorr = angle(sum(xcorr));

deltaF = phase_xcorr./2.0./pi*15000;

freqoffset_cp = deltaF;
        
