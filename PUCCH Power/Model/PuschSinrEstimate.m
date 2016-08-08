

function [Sinr] = PuschSinrEstimate(estChannelGrid,drsindex,noiseEst)

q = size(estChannelGrid,3);

Hp =  estChannelGrid(drsindex,[4 11],:);

Sp = mean(mean(mean(Hp.*conj(Hp))));

Sinr = 10*log10(Sp/noiseEst);




