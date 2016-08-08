
clc;
clear;


evm = 0.005:0.005:0.50
sinr = 1./(evm.*evm);
sinr_db = 10*log10(sinr);

figure();
plot(evm,sinr_db,'r-+');
grid on;
xlabel('EVM');
ylabel('SINR(dB)')
title('EVM to SINR');



