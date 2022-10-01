% First Computer project
% Student name: Reza SoleimaniMajd
% Student ID: 100880218
clear;
clc;
%% Task one

% Generate a set of 1000 i.i.d. random numbers
X = 4 * rand(1000, 1) - 2;

% Display histogram of the data:
figure(1);
h = histogram(X);
title('Uniform distribution of 1000 data samples')
xlabel('Bins')
ylabel('value')

% Calculate the pdf:
pdf = h.Values/1000;

% Calculate the CDF:
CDF = cumsum(pdf);

% Plot CDF:
figure(2)
plot(CDF)

% Another way of plotting:
% Since this way is more accurate, we plot this as the final result:
sorted = sort(X);
plot(sorted)
title('CDF of uniform r.v.')
xlabel('x')
ylabel('FX(x)')

%% Task two

% Generate a set of 1000 i.i.d random Gaussian numbers
G_x = randn(1000, 1);

% Display the histogram of data:
figure(3);
G_h = histogram(G_x);
title('Gaussian distribution of 1000 data samples')
xlabel('Bins')
ylabel('value')

% Calculate the pdf:
G_pdf = G_h.Values/1000;

% Calculate the CDF:
G_CDF = cumsum(G_pdf);

% Plot CDF:
figure(4);
plot(G_CDF)
title('CDF of Gaussian r.v.')
xlabel('x')
ylabel('FX(x)')

%% Task three

% Calculate autocorrolation funtion

% We can also use autocorr() function built in MATLAB 
% acf = autocorr(X, NumLags=50);
acf = Rx_est(X, 50);
figure(5)
stem(acf, "filled");
title('Autocorrolation of the random process')
xlabel('M (number of lags, or taw)')
ylabel('R(taw)')

% Calculate power density spectrum
S_x = fft(acf);
figure(6);
% Since we have 51 datapoints
fshift = -25:25;
stem(fshift, abs(fftshift(S_x)), "filled");
title('Power density of Auto-corrrolation function')
xlabel('f (Hz)')
ylabel('Power (watt)')

%% Task four

% Calculate the output of the linear filter:
Y = zeros(size(X));
Y_n_1 = 0;
for i=1:length(X)
    Y(i) = 0.95 * Y_n_1  + X(i);
    Y_n_1 = Y(i);
end

% Calculate autocorrolation funtion of filtered random process

% acf_filtered = autocorr(Y, NumLags=50);
acf_filtered = Rx_est(Y, 50);
figure(7)
stem(acf_filtered, "filled");
title('Autocorrolation of the filtered random process')
xlabel('M (number of lags, or taw)')
ylabel('R(taw)')

% Calculate power density spectrum of filtered random process
S_y = fft(acf_filtered);
figure(8);
fshift = -25:25;
stem(fshift, abs(fftshift(S_y)), "filled");
title('Power density of filtered auto-corrrolation function')
xlabel('f (Hz)')
ylabel('Power (watt)')

%% Task five

% Generate X sequence
Binary_X = ones(10^8,1);

% Randomly select half of them. Since P(x=1) = P(x=-1) = 0.5
% and then, change those from 1 to -1
random_idx = randperm(10^8);
idx = random_idx(1:10^8/2);
Binary_X(idx) = -1;

% Genarate N
N = randn(10^8,1);

% Generate Y
Binary_Y = Binary_X + N;

% Estimation (using sign function)
X_hat = sign(Binary_Y);

% Count the error:

% Since the error is abs(1--1), so it should be divided by 2 to count the 
% number of error occurance 
Error = sum(abs(X_hat - Binary_X))/2;
Error_per = 100 * Error/10^8;
disp("Error Percentage is:")
disp(Error_per)


%% Increse SNR
SNR = 0:1:10;
sigma = sqrt(1./10.^(SNR/10));

Error_per_SNR = zeros(size(sigma));
for i=1:numel(sigma)
    % Generate
    Noise = normrnd(0,sigma(i), [10^8,1]);
    Noisy_output = Binary_X + Noise;
    
    % Estimate
    X_hat = sign(Noisy_output);

    % Count the error:
    Error = sum(abs(X_hat - Binary_X))/2;
    Error_per_SNR(i) = 100 * Error/10^8;
end

figure(9)
semilogy(SNR, Error_per_SNR)
grid on
title('Error percentage for different SNR')
xlabel('SNR (dB)')
ylabel('Error')
