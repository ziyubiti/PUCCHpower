


% alpha low pass filter ver
clear;
clc;

alpha = [ 0.1 0.25 0.5 0.75 0.9];

data = random('Normal',0,100,[1,100]);

for i = 1:size(alpha,2)
    coee = alpha(i);
    for j = 1:size(data,2)
        if j == 1
            D(i,j) = data(1);
        else
            D(i,j) = coee*D(i,j-1) + (1-coee)*data(j);
        end;
    end;
end;

figure();
plot(data,'-bo');
grid on;
hold on;
plot(D(1,:),'-r*');
hold on;
plot(D(2,:),'-rp');
hold on;
plot(D(3,:),'-c>');
hold on;
plot(D(4,:),'-c<');
hold on;
plot(D(5,:),'-g^');
