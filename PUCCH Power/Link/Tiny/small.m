

% small distance los pathloss

clear;

%% Alt 1: free space
F = 2300;   % MHz;
D = 1:100:36000;  % Km

PL = 20*log10(F) + 20*log10(D) + 32.45;   % free space,

figure();
plot(D,PL,'r-p');
grid on;

%% Alt 2:cost-231 Hata model
f1 = 600;
f2 = 2300;
d = 0.05:0.01:10.0;
ht = 2.0;
hr = 2.0;
a_hr = 0
Cm = -23.6;    % rural
PL1 = 46.3+33.9*log10(f1)-13.82*log10(ht)-a_hr+(44.9-6.55*log10(ht))*log10(d)+Cm;
PL2 = 46.3+33.9*log10(f2)-13.82*log10(ht)-a_hr+(44.9-6.55*log10(ht))*log10(d)+Cm;
% pl = 46.33+(44.9-6.55*log10(ht))*log10(d)+33.9*log10(f)-((1.1*log10(f)-0.7)*hr-1.56*log10(f)+0.8)-13.82*log(ht)
pl1 = 46.33+(44.9-6.55*log10(ht))*log10(d) + 33.9*log10(f1)-((1.1*log10(f1)-0.7)*hr-1.56*log10(f1)+0.8)-13.82*log(ht);
pl2 = 46.33+(44.9-6.55*log10(ht))*log10(d) + 33.9*log10(f2)-((1.1*log10(f2)-0.7)*hr-1.56*log10(f2)+0.8)-13.82*log(ht);

pl_los = 42.64+20*log10(f2)+26*log10(d);
% R = 6378;
% PL2_earth = 46.3+33.9*log10(f2)-13.82*log10(ht+d*sin(asin(d/R)))-a_hr+(44.9-6.55*log10(ht))*log10(d)+Cm;


figure();
% plot(d*1000,PL1,'r-',d*1000,PL2,'b-');
plot(d*1000,PL1,'b-',d*1000,pl_los,'r-');
grid on;
title('Ht = 2m');
xlabel('Distance(m)');
ylabel('PathLoss(dB)');

figure();
plot(d,pl_los,'r-o');
grid on;



figure();
surf(d,ht,PL2);


figure();
plot(d*1000,pl1,'r-',d*1000,PL1,'b-');
grid on;
title('Ht = 2m');
xlabel('Distance(m)');
ylabel('PathLoss(dB)');





%% Alt 3:sea far distance
f = 1800;               % MHz
d = 10:0.1:50;        %km
ht = 100.0;
hr = 4;
a = 5;
Lboat = 0;
Learth = 0;

lamda = 300/f;
L0 = 20*log10(f) + 20*log10(d) + 32.45;   % free space,
PL = L0 - 10*log10(2-2*cos(4*pi*ht*hr./(1000*d*lamda))) + a + Lboat + Learth;



figure();
% plot(d,PL,'r-',d,L0,'b-');
plot(d,PL,'r-');
grid on;
title('Ht = 100m');
xlabel('Distance(km)');
ylabel('PathLoss(dB)');


%%







