


% pusch xcorr ver

% between board IQ data and ideal Tx data;

clear;
clc;

% board data;
dir = 'I:\20150709\2';

File_name = strcat(dir,'\TD_UL_LowHigh_40ms.txt');

A = load(File_name);

B = reshape(A,2,[]);

C = B.';

D = C(:,1) + 1j*C(:,2);

E = abs(D);


N = 2;  % 0 1 2 3
figure();
plot(E(1+307200*N:307200+307200*N,1));
grid on;

F = D(1+307200*N:307200+307200*N,1)/2.^15;            % 10 ms TD IQ  board data;

File_name2 = strcat(dir,'\FD_Rx.mat');
G = load(File_name2);
FD_rx = G.rxSubframe;

% sf2 ideal data;
A = load('A1.3_FDgrid.mat');
B = load('A1.3_TD.mat');

FD_data = A.txSubframe;
TD_data = B.txWaveform;

% figure();
% plot(real(FD_data(:,4)));
% grid on;
% figure();
% plot(imag(FD_data(:,4)));
% grid on;
% 
% 
% figure();
% plot(abs(FD_data(:,4)));
% grid on;

% dmrs xcorr
DMRS_pos = 2048*3+160+144*2;
DMRS_TD = TD_data(DMRS_pos+1:DMRS_pos+2048+144,1);

DMRS_xcorr = circle_corr(F.',DMRS_TD.');

figure();
plot(abs(DMRS_xcorr));
grid on;

% sf2 xcorr
sf2_xcorr = circle_corr(F.',TD_data.');

figure();
plot(abs(sf2_xcorr));
grid on;


% FD xcorr, 1st symbol
sym1_xcorr = circle_corr(FD_rx(:,1).',FD_data(:,1).');

figure();
plot(abs(sym1_xcorr));
grid on;


