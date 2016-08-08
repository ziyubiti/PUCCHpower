
% File_name:   *.bin
% flag: 'LLHH'  ‘HHLL’
% baijie 2015.4



function [output] = File_bin(File_name,flag)

dir = 'D:\2016 ver\0808 128';


File_name = strcat(dir,'\gRadIqBuf_Dl_Ant0.bin');

flag = 'LLHH';  % HHLL, LLHH;       intel board is LLHH;  QC waveform also
% flag = 'HHLL'; 
Port2_flag = 0;   % IQ数据是单天线还是串行天线0/1

fileID = fopen(File_name);
A = fread(fileID);                % byte 
fclose(fileID);

B = (reshape(A,2,[])).';           % 16 bit,  2 byte for 2 column;

if (strcmp(flag,'LLHH') == 1)
    C = B(:,2).*2.^8 + B(:,1);
elseif (strcmp(flag,'HHLL') == 1)
    C = B(:,1).*2.^8 + B(:,2);
end;

D = C(1:2:end,1) + 1j*C(2:2:end,1);
E = TypeConv(D);






Ifile = strcat(dir, '\TD_LowHigh.txt');

fid1 = fopen(Ifile,'w+');
if (Port2_flag == 1)
    for i = 1:1:0.5*size(E,1)
        fprintf(fid1,'%d  \n%d \n' ,real(E(i,1))+real(E(i+0.5*size(E,1),1)),imag(E(i+0.5*size(E,1),1))+imag(E(i,1)));
    end;
else    
    for i = 1:1:size(E,1)
        fprintf(fid1,'%d  \n%d \n' ,real(E(i,1)),imag(E(i,1)));
    end;
end;

fclose(fid1);
% 



figure();
plot(abs(E));
grid on;

figure();
plot(abs(E(1:1228800,1)));
grid on;


SFNum = 40;
Ifile2 = strcat(dir, '\TD_LowHigh_40ms.txt');

fid2 = fopen(Ifile2,'w+');
for i = 1:30720*SFNum
    fprintf(fid2,'%d  \n%d \n' ,real(E(i,1)),imag(E(i,1)));    
end;
fclose(fid2);



%%
function output = TypeConv(input)

Seq_Length = length(input);
output = zeros(Seq_Length,1);

for ii = 1:Seq_Length,
    if (real(input(ii)) < 32768),
        RePart = real(input(ii));
    else
        RePart = -(32767 - ((real(input(ii))-32768)-1));    % -(2.^16 - x);
    end;
    
    if (imag(input(ii)) < 32768),
        ImPart = imag(input(ii));
    else
        ImPart = -(32767 - ((imag(input(ii))-32768)-1));
    end;
    output(ii) = RePart + 1j*ImPart;
    %         output(ii) = 1j*RePart + ImPart;
end;

    
    



