
% dsp data only 1 frame repetition

clear;


A = load('D:\2015 work\01 RAN4 Test\N7625B testvector\TD_LowHigh.txt');

B = reshape(A,2,[]);

C = B.';

D = C(:,1) + 1j*C(:,2);

%wangxiaoyu
% start = 184320;
%D = D(start+1:1:start+307200*2);



E = abs(D);
sfnum = size(E,1)/30720;
% 
% figure();
% plot(E);
% grid on;
% 
% 
% 
figure();
plot(E(1:307200*4,1));
%plot(E);
grid on;

figure();
semilogy(E(1:307200*4,1));
grid on;

N = 0;
figure();
plot(E(1+307200*N:307200+307200*N,1));
grid on;
figure();
semilogy(E(1+307200*N:307200+307200*N,1));
grid on;


D2 = D(1+307200*N:307200+307200*N,1);        % 1 frame,start with sf9;


TM = lteTestModel('1.1','20MHz',1,'TDD');
 TM.SSC = 7;
 TM.TDDConfig = 1;

[tmwaveform,tmgrid,tmconfig] = lteTestModelTool(TM);   

%生成波形文件，二进制，用于高通平台
tt = round(tmwaveform*(2.^15));
tt2 = [real(tt) imag(tt)];
tt3 = reshape(tt2.',[],1);

Ifile = strcat(dir,'\ETM11B.bin');

fid1 = fopen(Ifile,'wb');
COUNT = fwrite(fid1,tt3,'int16');
fclose(fid1);



figure();
scatter(real(reshape(tmgrid,[],1)),imag(reshape(tmgrid,[],1)));
grid on;


figure();
scatter(real(tmgrid(:,17)),imag(tmgrid(:,17)));
grid on;

% pdsch_index = ltePDSCHIndices(tmconfig);

% unit test,prb positon; etm3.2
tmgrid = rxGrid;
PRBQ = [];
fd_data = tmgrid(:,[7:14:end]);
for iSF = 1:10
    if (iSF ~= 6 & iSF ~= 4 & iSF ~= 5)         % dsp data, start sf 9
        tmp_data = fd_data(:,iSF);
        tmp2 = reshape(tmp_data,12,[]);
        tmp3 = tmp2.*conj(tmp2);
        tmp4 = mean(tmp3,1);
        tmp5 = 10*log10(tmp4);
        tmp6 = sign(tmp5);
        prb_qam = find(tmp6==-1)-1;
        PRBQ = [PRBQ prb_qam.'];
    end;
end;

tmp_td = tmwaveform;

freqoffset_cp = 2000;            % unit test;
% big freq offset compensation
nts = 0:size(tmp_td,1)-1;
Ts = 1e-6./30.72;
comp = exp(-j*2*pi*freqoffset_cp*nts*Ts);
tmp_td = tmp_td.*comp.';

rxSubframe = lteSCFDMADemodulateV2(frc,tmp_td);

% tmwaveform = waveform;
D2 = tmp_td;

rxGrid = lteOFDMDemodulate(tmconfig,D2);

rxGrid = lteOFDMDemodulate(tmconfig,D2)/sqrt(2048)/2.^13*sqrt(2);

figure();
scatter(real(rxGrid(:,31-14)),imag(rxGrid(:,31-14)));
grid on;


figure();
scatter(real(rxGrid(:,7)),imag(rxGrid(:,7)),'bo');
grid on;
hold on;
scatter(real(rxGrid([2:6:end],1)),imag(rxGrid([2:6:end],1)),'r*');
hold on;
scatter(real(rxGrid([570:630],31-14)),imag(rxGrid([570:630],31-14)),'cp');  % PSS
hold on;
scatter(real(rxGrid([570:630],28-14)),imag(rxGrid([570:630],28-14)),'g>');  % SSS
hold on;
scatter(real(rxGrid([1:6:end],4)),imag(rxGrid([1:6:end],4)),'bo');    % control


%% test vector 
% rxGrid = [tmgrid(:,[127:140]) tmgrid(:,[1:126])];  % frame 0, SF9/SF0-SF8;
% rxGrid = [tmgrid(:,[267:280])  tmgrid(:,[141:266])];  % frame 1, SF19/SF10-SF18

%%
figure();
plot(abs(reshape(rxGrid,[],1)));
grid on;

% find the prb position
for i = 0:1:9
    symbol_idx = 1 + i*14:(i+1)*14;
    fd_sf = rxGrid(:,symbol_idx);
    figure(i+1);
    plot(abs(fd_sf));
    grid on;
end;

% nPRB index
% nPRB = [51 71 99 0 0 0 85 35 4 8 ];          % sf 9 /0 、1，5-8, to fit in dsp Data;
nPRB = [51 71 99 0 0 0 85 35 4 8 ];          % frame 0;
% nPRB = [76 75 93 0 0 0 34 20 53 23];   % frame 1;SF19/SF10-SF18

% obtain pdsch prb
fd_sf = [];
tmp_sf = [];
crs_idx = lteCellRSIndices(tmconfig);
for i = 0:1:9
     sf_idx = mod(i+9,10);
    symbol_idx = 1 + i*14:(i+1)*14;
    rxGrid2 = reshape(rxGrid(:,symbol_idx),[],1);
    rxGrid2(crs_idx) = 99;
    rxGrid3 = reshape(rxGrid2,1200,[]);
    sc_idx = [nPRB(i+1)*12+1: (nPRB(i+1)+1)*12];
    tmp_sf = rxGrid3(sc_idx,:);
    tmp_sf(:,1) = 0;
    if sf_idx == 1
        tmp_sf(:,[12 13 14]) = 0;
    end;
    if (sf_idx == 2 | sf_idx == 3 |sf_idx == 4)
        tmp_sf(:,:) = 0;
    end;
    fd_sf = [fd_sf tmp_sf];    
end;

Pm_sf = abs(fd_sf).*abs(fd_sf);

figure();
plot(reshape(Pm_sf,[],1));
grid on;

for i = 0:1:9
     sf_idx = mod(i+9,10);
    symbol_idx = 1 + i*14:(i+1)*14;
    Pm = reshape(Pm_sf(:,symbol_idx),[],1);
    pdsch_idx = find(Pm > 0 & Pm <50);
    avePm(i+1) = mean(Pm(pdsch_idx));
    avePm_dB(i+1) = 10*log10(avePm(i+1));
end;
    
figure();
stem(avePm_dB);
grid on;
title('SF9 SF0-SF8');
ylabel('PO(dB)');

close all;

% pcfich_index = ltePCFICHIndices(tmconfig);
% phich_index = ltePHICHIndices(tmconfig,{'ind','re'});