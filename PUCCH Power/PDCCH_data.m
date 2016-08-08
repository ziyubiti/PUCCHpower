clear;

A= load('E:\technote20150415\mipan\PDSCH - 1c1u\Debug\outputfiles\time_data9_interferenceData_PCI0.txt');

B = A(:,1) + 1j*A(:,2);

figure();
plot(abs(B));
grid on;

% cat1:  same CRS;      
% cat2:  diff CRS;
% cat3:  same CRS + data;
% cat4:  diff CRS + data

% pri: no 0 5 -10 10 ;


%% 100sf ver
snr = [   -6        -4        -2         0         2         4         6           8        10 ];

% 1c1u
bler1=[ 0.38      0.05         0         0         0 0 0 0 0];

% inter: -10dB, same CRS
 bler2=[   0.41      0.07         0         0         0 0 0 0 0 ];

% inter: 0dB, same CRS
  bler3=[     0.54      0.15         0         0         0  0 0 0 0];

% inter: 4.77dB,same CRS
bler4 = [   0.64      0.32      0.05      0.02         0         0         0           0         0     ];

% inter: 10dB,same CRS
bler5 = [      0.76      0.63      0.51      0.38       0.2      0.12      0.07       0.06      0.04  ];


figure(1);
semilogy(snr,bler1,'r-*',snr,bler2,'b-o',snr,bler3,'c-<',snr,bler4,'b-p',snr,bler5,'r->');
grid on;
title('PDCCH performace');
legend('No Interference','interference: -10dB, same CRS','interference: 0dB, same CRS','interference: 5dB, same CRS','interference: 10dB, same CRS',2);
xlabel('SNR');
ylaber('BLER');

%%

% inter: -10dB, diff CRS


% inter: 0dB, diff CRS


% inter: 10dB,diff CRS