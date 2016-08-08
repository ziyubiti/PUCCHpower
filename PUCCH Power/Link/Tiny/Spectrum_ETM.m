
clear;
  
tm = '3.1';     % Test model number
bw = '20MHz';     % Bandwidth

[timeDomainSig, txGrid, txInfo] = lteTestModelTool(tm,bw);

figure(1);
plot(abs(timeDomainSig));
grid on;


figure(2);
plot(abs(reshape(txGrid,[],1)));
grid on;

figure(3);
scatter(real(reshape(txGrid,[],1)),imag(reshape(txGrid,[],1)));
grid on;

% hPlotDLResourceGrid(txInfo,txGrid);

%% Plot Spectrogram
% Plot a spectrogram of the time domain signal.

% Compute spectrogram
[y,f,t,p] = spectrogram(TV_data, 2048, 0, 2048, txInfo.SamplingRate);  
[y,f,t,p] = spectrogram(timeDomainSig, 2048, 0, 2048, txInfo.SamplingRate);   %bj:中间重叠1024点，分段长2048；默认值不重叠

% Re-arrange frequency axis and spectrogram to put zero frequency in the
% middle of the axis i.e. represent as a complex baseband waveform
f = (f-txInfo.SamplingRate/2)/1e6;
p = fftshift(10*log10(abs(p)));

% Plot spectrogram，3D
figure;
surf(t*1000,f,p,'EdgeColor','none');   
xlabel('Time (ms)');
ylabel('Frequency (MHz)');
zlabel('Power (dB)');
title('Spectrogram of Test Model E-TM1.1, 20MHz');

%% bj: 2D plot
figure();
plot(f,p(:,1));grid on;


figure();
plot(10*log10(fftshift(abs(y(:,1)))));grid on;


figure();
tmp = fftshift(10*log10(abs(fft(timeDomainSig(1:2048,1)))));
plot(f,tmp);grid on;

h = hamming(2048);
figure();
tmp = fftshift(10*log10(abs(fft(timeDomainSig(1:2048,1)).*fft(h))));
plot(f,tmp);grid on;
