
% CRS channel estimation to calc time and freq offset
% s2 = s1 * exp(j*2pi*deltaF*delataT*nTs);
% normalized timeoffset ,express Ts;
% absolute freqoffset,express Hz;   
% small freq offset, [-1000, 1000];
% hest:4D, 1200*14*Nrx*Ntx

function [timeoffset,freqoffset] = lteDlTimeFreqOffsetEstimate(hest)

H = [];
H07 =  hest(:,[1 8],1,1);
H411 =  hest(:,[5 12],1,1);

H = [H07;H411];

Xcorr = conj(H(:,1)).*H(:,2);

phase_xcorr = angle(sum(Xcorr));

deltaF = phase_xcorr./2.0./pi/7.5*15000;
freqoffset = deltaF;
% 
% Ht = [];
% for n = 1:q
%     Ht = [Ht Hp(:,:,n)];
% end;
% 
% m = 6;
% 
% Xcorr2 = Ht([1:end-m],:).*conj(Ht([m+1:end],:));
% 
% phase_xcorr2 = angle(sum(sum(Xcorr2)));
% 
% deltaT = phase_xcorr2./2.0./pi/m*2048;
%  
% timeoffset = deltaT;

