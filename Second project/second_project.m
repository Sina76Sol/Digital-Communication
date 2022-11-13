% Second Project
% Student name: Reza Soleimani Majd
% Student ID: 100880218
clear;
clc;
%% ::::::::Task 1.1: Building the constellation for QPSK::::::::
c_QPSK = [sqrt(2)/2+sqrt(2)/2i, -sqrt(2)/2+sqrt(2)/2i,...
         -sqrt(2)/2-sqrt(2)/2i, sqrt(2)/2-sqrt(2)/2i];
QPSK_sigpower = pow2db(mean(abs(c_QPSK).^2));
M = length(c_QPSK);

% create original and modulated symbol
QPSK_x = randi([0 M-1],10^8,1);
QPSK_mod_x = genqammod(QPSK_x,c_QPSK);

% create awgn received symbol
QPSK_y = awgn(QPSK_mod_x,12,QPSK_sigpower);

%% Plot
h = scatterplot(QPSK_y);
hold on
scatterplot(c_QPSK,[],[],'r*',h)
grid

%% Experimental error rate for QPSK
SNR = 0:3:15;
QPSK_error_rate = zeros(size(SNR));

for i = 1:length(SNR)
    noisy_sig = awgn(QPSK_mod_x,SNR(i),QPSK_sigpower);
    QPSK_detected = genqamdemod(noisy_sig,c_QPSK);
    [numErrors,ser] = symerr(QPSK_x,QPSK_detected);
    QPSK_error_rate(i) = ser;
end

%% plot semilog graph
close;
semilogy(SNR, QPSK_error_rate,'-o');
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

%% ::::::::Tasl 1.2: Build the constellation for 4PAM::::::::
delta  = sqrt(4/5);
c_4pam = [-3*delta/2, -delta/2, delta/2, 3*delta/2];
PAM_sigpower = pow2db(mean(abs(c_4pam).^2));
M = length(c_4pam);

% create original and modulated symbol
PAM_x = randi([0 M-1],10^8,1);
PAM_mod_x = genqammod(PAM_x,c_4pam);

% create awgn received symbol
PAM_y = awgn(PAM_mod_x,12,PAM_sigpower);

%% Plot 4PAM
close;
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

%% plot semilog graph
close;
semilogy(SNR, PAM_error_rate,'-o')
hold on

%% Theoritical error rate for 4PAM
delta =sqrt(4/5);
Bit_energy = 0.5;
sigma = sqrt(Bit_energy./(10.^(SNR./10)));
PAM_Pe = 1.5*qfunc(delta./(2*sigma));
semilogy(SNR, PAM_Pe, 'r-*')
xlabel('SNR')
ylabel('error rate')
legend("experimental", "theoritical")
hold off

%% plot both graphs together
close
semilogy(SNR, PAM_error_rate)
hold on
semilogy(SNR, QPSK_error_rate)
xlabel('SNR')
ylabel('error rate')
legend("4PAM error rate", "QPSK error rate")

%% ::::::::Task 2.1: build the designed 4M constellation::::::::

% assign the constellation points
c_4design = [-sqrt(1.1), sqrt(0.9)*1i, -sqrt(0.9)*1i, sqrt(1.1)];
M4_sigpower = pow2db(mean(abs(c_4design).^2));
M = length(c_4design);

% create original and modulated symbols
M4_x = randi([0 M-1],10^4,1);
M4_mod_x = genqammod(M4_x,c_4design);

% create received symbols
M4_y = awgn(M4_mod_x,15,M4_sigpower);

%% plot the 4M designed  constellation
close;
h = scatterplot(M4_y);
hold on
scatterplot(c_4design,[],[],'r*',h)
grid;
hold off

%% calculating the theoritical upper bound for 4M
close
SNR = 0:3:15;
d_min = sqrt(2);
Bit_energy = db2pow(M4_sigpower)/log2(M);
sigma = sqrt(Bit_energy./(10.^(SNR./10)));
theory_ub = (M-1)*qfunc(d_min./(2*sigma));

%% simulation error rate for 4M design
M4_error_rate = zeros(size(SNR));
for i = 1:length(SNR)
    M4_noisy_sig = awgn(M4_mod_x,SNR(i),M4_sigpower);
    M4_detected = genqamdemod(M4_noisy_sig,c_4design);
    [M4_numErrors,M4_ser] = symerr(M4_x,M4_detected);
    M4_error_rate(i) = M4_ser;
end

%% plot the error rate for 4M designed:

close
semilogy(SNR, M4_error_rate,'-o')
hold on
semilogy(SNR, theory_ub,'-*')
xlabel('SNR')
ylabel('error rate')
legend("simulation error", "theory upper bound error")

%% ::::::::Task 2.2 build the designed 8M constellation::::::::
% assign the constellation points
c_8QAM = [0.5+0.5*1i, 0.5-0.5*1i, -0.5+0.5*1i, -0.5-0.5*1i...
    sqrt(3/2), -sqrt(3/2), sqrt(3/2)*1i, -sqrt(3/2)*1i];
QAM_sigpower = pow2db(mean(abs(c_8QAM).^2));
M = length(c_8QAM);

% create original an modulated symbol
x_8QAM = randi([0 M-1],10^4,1);
mod_x_8QAM = genqammod(x_8QAM,c_8QAM);

% create received symbols
QAM_y = awgn(mod_x_8QAM,12,QAM_sigpower);
%% Plot 8QAM constellation
close;
h = scatterplot(QAM_y);
hold on
scatterplot(c_8QAM,[],[],'r*',h)
grid;
hold off
%% Calculate theoritical upper bound
d_min = sqrt((sqrt(3/2)-0.5)^2+0.5^2);
Bit_energy = db2pow(QAM_sigpower)/log2(M);
sigma = sqrt(Bit_energy./(10.^(SNR./10)));
theory_ub_8QAM = (M-1)*qfunc(d_min./(2*sigma));

%% Simulation error for 8QAM
QAM_error_rate = zeros(size(SNR));
for i = 1:length(SNR)
    QAM_noisy_sig = awgn(mod_x_8QAM,SNR(i),QAM_sigpower);
    QAM_detected = genqamdemod(QAM_noisy_sig,c_8QAM);
    [QAM_numErrors,QAM_ser] = symerr(x_8QAM,QAM_detected);
    QAM_error_rate(i) = QAM_ser;
end
%% plot error rate for 8QAM design
close
semilogy(SNR, QAM_error_rate,'-o')
hold on
semilogy(SNR, theory_ub_8QAM,'-*')
xlabel('SNR')
ylabel('error rate')
legend("simulation error", "theory upper bound error")

%% ::::::::design 8PKS constellation:::::::
% assign the constellation points
c_8PSK = [1, -1, 1i, -1i, ...
          sqrt(1/2)+sqrt(1/2)*1i, -sqrt(1/2)+sqrt(1/2)*1i,...
          -sqrt(1/2)-sqrt(1/2)*1i, sqrt(1/2)-sqrt(1/2)*1i];
OPSK_sigpower = pow2db(mean(abs(c_8PSK).^2));
M = length(c_8PSK);

% create original an modulated symbol
x_8PSK = randi([0 M-1],10^4,1);
mod_x_8PSK = genqammod(x_8PSK,c_8PSK);

% create received symbols
OPSK_y = awgn(mod_x_8PSK,12,OPSK_sigpower);
%% plot 8PSK constellation
close;
h = scatterplot(OPSK_y);
hold on
scatterplot(c_8PSK,[],[],'r*',h)
grid;
hold off

%% Simulation error for 8PSK
OPSK_error_rate = zeros(size(SNR));
for i = 1:length(SNR)
    OPSK_noisy_sig = awgn(mod_x_8PSK,SNR(i),OPSK_sigpower);
    OPSK_detected = genqamdemod(OPSK_noisy_sig,c_8PSK);
    [OPSK_numErrors,OPSK_ser] = symerr(x_8PSK,OPSK_detected);
    OPSK_error_rate(i) = OPSK_ser;
end
%% compare 8PSK and 8QAM error rate
close
semilogy(SNR, QAM_error_rate,'-o')
hold on
semilogy(SNR, OPSK_error_rate,'-*')
xlabel('SNR')
ylabel('error rate')
legend("8QAM error", "8PSK error")
