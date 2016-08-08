% binary file for QC waveform generation

clear;

% LTE
TM = lteTestModel('1.1','20MHz',1,'TDD');
[tmwaveform,tmgrid,tmconfig] = lteTestModelTool(TM);   

m1 = max(real(tmwaveform));
m2 = max(imag(tmwaveform));


% single tone
f = 17.1e6;         % cf
nts = 0:614400;     % 20ms
Ts = 1e-6./30.72;
tdata = exp(j*2*pi*f*nts*Ts);


%% filter to improve ACLR
fir = dsp.LowpassFilter();
fir.SampleRate = tmconfig.SamplingRate;
fir.PassbandFrequency = tmconfig.NDLRB*180e3/2;
fir.StopbandFrequency = 10e6;

% Apply filter
waveform = step(fir,tmwaveform*10);

m3 = max(real(waveform));
m4 = max(imag(waveform));

%% verify evm after filter
offset = lteDLFrameOffset(tmconfig,waveform,'TestEVM');
cec.PilotAverage = 'TestEVM';
evmsettings.EnablePlotting = 'Off';
evm = hPDSCHEVM(tmconfig,cec,waveform(1+offset:end,:),evmsettings);


%% 生成波形文件，二进制，用于高通平台
tt = round(waveform*(2.^15));
tt2 = [real(tt) imag(tt)];
tt3 = reshape(tt2.',[],1);

binfile = 'E:\technote20150415\PUCCH Power\Link\Tiny\LTE_ETM11_QC.bin';

fid1 = fopen(binfile,'wb');
COUNT = fwrite(fid1,tt3,'int16');
fclose(fid1);

matfile = 'E:\technote20150415\PUCCH Power\Link\Tiny\LTE_ETM11_QC_mat.bin';
fid2 = fopen(matfile,'w+');
for i = 1:1:size(waveform,1)
    fprintf(fid2,'%f       %f \n' ,real(waveform(i,1)),imag(waveform(i,1)));    
end;
fclose(fid2);

