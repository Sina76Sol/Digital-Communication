clear;
clc;
%% Building the constellation for QPSK
c_QPSK = [sqrt(2)/2+sqrt(2)/2i, -sqrt(2)/2+sqrt(2)/2i,...
         -sqrt(2)/2-sqrt(2)/2i, sqrt(2)/2-sqrt(2)/2i];
QPSK_sigpower = pow2db(mean(abs(c_QPSK).^2));
M = length(c_QPSK);

QPSK_x = randi([0 M-1],10^4,1);
QPSK_mod_x = genqammod(QPSK_x,c_QPSK);

QPSK_y = awgn(QPSK_mod_x,3,QPSK_sigpower);

h = scatterplot(QPSK_y);
hold on
scatterplot(c_QPSK,[],[],'r*',h)
grid

%% Experimental error rate
SNR = 0:3:15;
error_rate = zeros(size(SNR));
for i = 1:length(SNR)
    noisy_sig = awgn(QPSK_mod_x,SNR(i),QPSK_sigpower);
    QPSK_detected = genqamdemod(noisy_sig,c_QPSK);
    [numErrors,ser] = symerr(QPSK_x,QPSK_detected);
    error_rate(i) = ser;
end

semilogy(error_rate,'-o')
hold on
%% Theoritical error rate
delta =sqrt(2);
Bit_energy = 0.5;
sigma = sqrt(Bit_energy./(10.^(SNR./10)));
Pe = 1-(1-qfunc(delta./(2*sigma))).^2;
semilogy(Pe, 'r-*')
legend("experimental", "theoritical")

%% Build the constellation for 4PAM
delta  = sqrt(4/5);
c_4pam = [-3*delta/2, -delta/2, delta/2, 3*delta/2];
sigpower = pow2db(mean(abs(c_4pam).^2));
M = length(c_4pam);

QPSK_x = randi([0 M-1],10^4,1);
modData = genqammod(data,c_4pam);

rxSig = awgn(modData,3,sigpower);
h = scatterplot(rxSig);
hold on
scatterplot(c_4pam,[],[],'r*',h)
grid
