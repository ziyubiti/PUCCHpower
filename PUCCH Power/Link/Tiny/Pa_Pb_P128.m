

% power offset calc for 128

% para
R = 0.2:0.1:0.9;
Nprb = 100;
Nre = 8;


% calc
rho_crs = 6*R;
rho_A = 100.0/Nprb;
rho_B = 1200*(1-R)/(Nprb*Nre);


Pa_dB = 10*log10(100.0./(Nprb*6*R));
Pb = 12*(1-R)./Nre;

deltaAD = round(40*log10(rho_crs));

deltaAD_hex = -dec2hex(deltaAD);

result = [R;rho_crs;repmat(rho_A,[1,size(R,2)]);rho_B;Pa_dB;Pb;deltaAD].';

figure();
surf(result);
grid on;




