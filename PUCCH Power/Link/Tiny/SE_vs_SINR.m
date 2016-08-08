

% Spectrum efficiency vs SINR, uplink , 10RB

sinr = [0.365 	1.259 	1.840 	2.774 	3.340 	4.095 	4.682 	5.480 	5.998 	6.584 	...
    7.294 	7.818 	8.483 	9.223 	9.800 	10.680 	11.350 	11.695 	12.690 	13.650 	14.320 	14.710 	15.480 	16.320 	17.210 	18.080 	18.550 	19.250 	20.320 ];

mcs = 0:28;



figure();
plot(mcs,sinr,'r-o');
grid on;
title('UL MCS vs SINR');
xlabel('MCS');
ylabel('SINR for BLER 10%');


% 10RB
TBS_QPSK = [256 344 424 568 696 872 1032 1224 1384 1544 1736]; 
TBS_16QAM = [1736 2024 2280 2536 2856 3112 3240 3624 4008 4264] ;
TBS_64QAM = [4264 4584 4968 5352 5736 5992 6200 7480];

Ndata = 12*12*10;
SE_QPSK = (TBS_QPSK+24)/Ndata;
SE_16QAM = (TBS_16QAM+24)/Ndata;
SE_64QAM(1:6) = (TBS_64QAM(1:6)+24)/Ndata;
SE_64QAM(7:8) = (TBS_64QAM(7:8)+24*3)/Ndata;
sinr_QPSK = [0.365 	1.259 	1.840 	2.774 	3.340 	4.095 	4.682 	5.480 	5.998 	6.584 	 7.294 ];
sinr_16QAM = [7.818 	8.483 	9.223 	9.800 	10.680 	11.350 	11.695 	12.690 	13.650 	14.320];
sinr_64QAM = [14.710 	15.480 	16.320 	17.210 	18.080 	18.550 	19.250 	20.320];


figure();
plot(SE_QPSK,sinr_QPSK,'r-o',SE_16QAM,sinr_16QAM,'b-o',SE_64QAM,sinr_64QAM,'g-o');
grid on;
title('UL SE vs SINR');
xlabel('Spectrum Efficiency');
ylabel('SINR for BLER 10%');
legend('QPSK','16QAM','64QAM','Location','SouthEast');
















