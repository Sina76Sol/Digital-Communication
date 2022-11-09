clear;
clc;
%% Building the constellation for QPSK
c_QPSK = [sqrt(2)/2+sqrt(2)/2i, -sqrt(2)/2+sqrt(2)/2i,...
         -sqrt(2)/2-sqrt(2)/2i, sqrt(2)/2-sqrt(2)/2i];
QPSK_sigpower = pow2db(mean(abs(c_QPSK).^2));
M = length(c_QPSK);

QPSK_x = randi([0 M-1],10^8,1);
QPSK_mod_x = genqammod(QPSK_x,c_QPSK);

QPSK_y = awgn(QPSK_mod_x,15,QPSK_sigpower);

h = scatterplot(QPSK_y);
hold on
scatterplot(c_QPSK,[],[],'r*',h)
grid

%% Experimental error rate for QPSK
SNR = 0:3:15;
error_rate = zeros(size(SNR));
for i = 1:length(SNR)
    noisy_sig = awgn(QPSK_mod_x,SNR(i),QPSK_sigpower);
    QPSK_detected = genqamdemod(noisy_sig,c_QPSK);
    [numErrors,ser] = symerr(QPSK_x,QPSK_detected);
    error_rate(i) = ser;
end

semilogy(SNR, error_rate,'-o')
hold on
%% Theoritical error rate for QPSK
delta =sqrt(2);
Bit_energy = 0.5;
sigma = sqrt(Bit_energy./(10.^(SNR./10)));
Pe = 1-(1-qfunc(delta./(2*sigma))).^2;
semilogy(SNR, Pe, 'r-*')
xlabel('SNR')
ylabel('error rate')
legend("experimental", "theoritical")

%% Build the constellation for 4PAM
delta  = sqrt(4/5);
c_4pam = [-3*delta/2, -delta/2, delta/2, 3*delta/2];
PAM_sigpower = pow2db(mean(abs(c_4pam).^2));
M = length(c_4pam);

PAM_x = randi([0 M-1],10^8,1);
PAM_mod_x = genqammod(PAM_x,c_4pam);

PAM_y = awgn(PAM_mod_x,12,PAM_sigpower);
h = scatterplot(PAM_y);
hold on
scatterplot(c_4pam,[],[],'r*',h)
grid
%% Experimental error rate for 4-PAM
SNR = 0:3:15;
PAM_error_rate = zeros(size(SNR));
for i = 1:length(SNR)
    PAM_noisy_sig = awgn(PAM_mod_x,SNR(i),PAM_sigpower);
    PAM_detected = genqamdemod(PAM_noisy_sig,c_4pam);
    [PAM_numErrors,PAM_ser] = symerr(PAM_x,PAM_detected);
    PAM_error_rate(i) = PAM_ser;
end

semilogy(SNR, PAM_error_rate,'-o')
hold on
%% Theoritical error rate for QPSK
delta =sqrt(4/5);
Bit_energy = 0.5;
sigma = sqrt(Bit_energy./(10.^(SNR./10)));
PAM_Pe = 1.5*qfunc(delta./(2*sigma));
semilogy(SNR, PAM_Pe, 'r-*')
xlabel('SNR')
ylabel('error rate')
legend("experimental", "theoritical")
