
clear;

% RSRP/RSRQ vs. PRB occupy ratio

RQ = -11.0:0.1:-6.1;
RQ_liner = 10.^(RQ/10);

R = 1./RQ_liner./1.024./8-0.5;

figure();
plot(RQ,R,'r-o');
grid on;
title('PRB occupy Ratio vs RSRP/RSRQ');
xlabel('RSRQ(dB)');
ylabel('PRB occupy Ratio');

dataset = [RQ.' R.'];



