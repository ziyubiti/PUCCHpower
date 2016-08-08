close all;
clear;
path(path,'../../Model');

debug_flag = 0;
dir = 'D:\2016 ver\0804 bigpowerevm';
File_name = strcat(dir,'\TD_LowHigh.txt');

A = load(File_name);
B = reshape(A,2,[]);
C = B.';
D = C(:,1) + 1j*C(:,2);

E = abs(D);
sfnum = size(E,1)/30720;
% 
figure();
plot(E);
grid on;
% 
figure();
semilogy(E);
grid on;
% 
figure();
plot(E(1:307200*4,1));
grid on;

figure();
semilogy(E(1:307200*4,1));
grid on;

F = 10*log10(E.*E)-90.3;   %dBFS

figure();
plot(F(1:307200*4,1));
grid on;


dl_offset = 1;
N = 1;  % 0 1 2 3
figure();
plot(E(dl_offset+1+307200*N:dl_offset+307200+307200*N,1));
grid on;
figure();
semilogy(E(dl_offset+1+307200*N:dl_offset+307200+307200*N,1));
grid on;

D2 = D(1+307200*N:307200+307200*N,1);   % ant 0
% D2 = D(1+307200*N+3686400:307200+307200*N+3686400,1);   % ant1
% D3 = D2(161:2208,1);
% D4 = fftshift(fft(D3))/sqrt(2048)/2.^15;
% Dtmp = fft(D3)/sqrt(2048);

figure();
plot(abs(D2));
grid on;

figure();
plot(real(D2),'r-');
hold on;
plot(imag(D2),'b-');
grid on;

% plot
tmconfig.NDLRB = 100;
tmconfig.CellRefP = 1;
tmconfig.NCellID = 1;
D2_align = D(1+30720:307200+30720,1);
figure();
plot(abs(D2_align));
grid on;

rxGrid = lteOFDMDemodulate(tmconfig,D2_align);    % start sf9

rxGrid = rxGrid/sqrt(2048)/2.^13/2;   % big power ,more /2.

CRSIndices0 = lteCellRSIndices(tmconfig,0);
CRSIndices1 = lteCellRSIndices(tmconfig,1)-16800;

[crs1] = lteExtractResources(CRSIndices1, rxGrid);


figure();
scatter(real(reshape(crs1,[],1)),imag(reshape(crs1,[],1)));
grid on;

for i = 0:1:9
    figure();
    scatter(real(reshape(rxGrid(:,[1+14*i:14*(i+1)]),[],1)),imag(reshape(rxGrid(:,[1+14*i:14*(i+1)]),[],1)));
    grid on;
end;

figure();
plot(abs(reshape(rxGrid,[],1)));
grid on;

figure();
scatter(real(rxGrid(:,31-14)),imag(rxGrid(:,31-14)));
grid on;


figure();
scatter(real(rxGrid(:,7)),imag(rxGrid(:,7)),'bo');
grid on;
hold on;
scatter(real(rxGrid([2:6:end],1)),imag(rxGrid([2:6:end],1)),'r*');
hold on;
scatter(real(rxGrid([5:6:end],1)),imag(rxGrid([5:6:end],1)),'r*');
hold on;
scatter(real(rxGrid([570:630],31-14)),imag(rxGrid([570:630],31-14)),'cp');  % PSS
hold on;
scatter(real(rxGrid([570:630],28-14)),imag(rxGrid([570:630],28-14)),'g>');  % SSS
hold on;
scatter(real(rxGrid([1:6:end],4)),imag(rxGrid([1:6:end],4)),'bo');    % control



% unit test for FDD, FFT for 1st Symbol

t_data = D(161:2208,1);
% t_data = D2(161:2208,1);
% t_data = D2(62140+1:62140+2208,1);
fd_data = fftshift(fft(t_data));
figure();
stem(abs(fd_data));
grid on;

if debug_flag == 1
    % add freqoffset
    freqoffset_cp = 0;            % unit test;
    % big freq offset compensation
    nts = 0:size(t_data,1)-1;
    Ts = 1e-6./30.72;
    comp = exp(j*2*pi*freqoffset_cp*nts*Ts);
    t_data2 = t_data.*comp.';
    
    % oversampling
    t_data3 = resample(t_data2,3,1);
    
    f = 17.1e6;
    nts = 0:2047;
    Ts = 1e-6./30.72;
    tdata4 = exp(j*2*pi*f*nts*Ts);
    fd_data = fftshift(fft(tdata4));
    % fd_data = (fft(tdata4));
    figure();
    stem(abs(fd_data));            % position 117, 1025 is DC,  (117-1025)*15k +30.72M = 17.1M;  because fftshift.
    grid on;
    
    [y,f,t,p] = spectrogram(tdata4, 2048, 0, 2048, 30.72e6);
    figure();
    stem(f,p(:,1));
    grid on;
    
end;
% 
% figure();
% plot(abs(D4));
% grid on;
% 
% D5 = D4([425:824 826:1625],1);
% 
% figure();
% semilogy(abs(D5));
% grid on;



%%
TA = 704;
adjust_factor = 704;%823;%-298; 14155  -5605
sf = 2;              % start sf0, to get sf2,
% offset = TA + 30720 + sf*30720 - TA ;   % start with sf9
offset =  sf*30720 - TA ;     % start with sf0
UL_data = D2(offset+1+adjust_factor:offset+30720+adjust_factor,1);


figure();
plot(abs(UL_data));
grid on;

figure();
plot(abs(fftshift(fft(UL_data(2208+1:2208+2048,1)))));
grid on;


% 
% P_td = mean(UL_data.*conj(UL_data));
% 
% P_td_dB = 10*log10(P_td);
% 
% AD9361_gain = 35;
% LNA_gain = 18;
% 
% Puu_td_dBFS = P_td_dB -10*log10(2.^30);
% Puu = P_td_dB+13-4.1-AD9361_gain-LNA_gain -10*log10(2.^30);


% freq 1st 25PRB power
frc.NULRB = 100;
frc.CyclicPrefixUL = 'Normal';
fd_data_all = lteSCFDMADemodulateV2(frc, UL_data);

FDgrid = reshape(fd_data_all,[],1);

figure();
plot(abs(FDgrid));
grid on;

pwr_fd = fd_data_all.*conj(fd_data_all);      

figure();
semilogy(pwr_fd(:,1));
grid on;

NPRB = [0 0];       % start and end

P_fd = sum(pwr_fd(NPRB(1)*12+1:NPRB(end)*12+12,:),1);
P_fd_N0 = mean(pwr_fd(NPRB(end)*12+13:end,:),1)*300;

P_fd_dB = 10*log10(P_fd/2048);
P_fd_dB_N0 = 10*log10(P_fd_N0/2048);

AD9361_gain = 47;
LNA_gain = 15;

Puu_fd_dBFS = P_fd_dB -10*log10(2.^30);
Puu_fd_dBFS_N0 = P_fd_dB_N0 -10*log10(2.^30);

Puu_fd = P_fd_dB+13-4.1-AD9361_gain-LNA_gain -10*log10(2.^30);
Puu_fd_N0 = P_fd_dB_N0+13-4.1-AD9361_gain-LNA_gain -10*log10(2.^30);

figure();
plot(Puu_fd_dBFS);
grid on;

figure();
plot(Puu_fd);
grid on;







% maxAmp = max([max(abs(real(UL_data)))  max(abs(imag(UL_data)))]);
% n = log2(maxAmp);
% 
% UL_data2 = UL_data./2.^(ceil(n));

UL_data2 = UL_data./2.^15;    % Q15


Ifile = strcat(dir,'\UL_SF2.txt');


fid1 = fopen(Ifile,'w+');
for i = 1:1:size(UL_data2,1)
    fprintf(fid1,'%f       %f \n' ,real(UL_data2(i,1)),imag(UL_data2(i,1)));    
end;
fclose(fid1);

close all;


%% IQ bit 
% signed 16 bit ,max 2.^15 = 32768;
% 
% n = 1:15;
% np = 4.^n;
% np_dBFS = 10*log10(np/2.^30);
% 
% figure();
% plot(np_dBFS,'r-o');
% grid on;






