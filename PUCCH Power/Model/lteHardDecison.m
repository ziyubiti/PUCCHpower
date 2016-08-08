function [HardBit] = lteHardDecison(rxEncodedBits)

len = length(rxEncodedBits);
HardBit = zeros(len,1);
for i = 1:len
    if (rxEncodedBits(i) < 0)
        HardBit(i) = 0;
    else
        HardBit(i) = 1;
    end;
end;
