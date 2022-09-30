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

% Get the value of each bin:
disp("Bin values are:");
disp(h.Values);

% Calculate the pdf:
pdf = h.Values/1000;

% Calculate the CDF:
CDF = cumsum(pdf);

% Plot CDF:
figure(2)
plot(CDF)

% Another way of plotting:
sorted = sort(X);
plot(sorted)

%% Task two

% Generate a set of 1000 i.i.d random Gausian numbers
G_x = randn(1000, 1);

% Display the histogram of data:
figure(3);
G_h = histogram(G_x);

% Get the value of each bin:
disp("Bin values are:");
disp(G_h.Values);

% Calculate the pdf:
G_pdf = G_h.Values/1000;

% Calculate the CDF:
G_CDF = cumsum(G_pdf);

% Plot CDF:
figure(4);
plot(G_CDF)

%% Task three

% Calculate autocorrolation funtion
acf = autocorr(X, NumLags=50);
figure(5)
stem(acf, "filled");

% Calculate power density spectrum
S_x = fft(acf);
figure(6);
fshift = -25:25;
stem(fshift, abs(fftshift(S_x)), "filled");

%% Task four

% Calculate the output of the linear filter:
Y = zeros(size(X));
Y_n_1 = 0;
for i=1:length(X)
    Y(i) = 0.95 * Y_n_1  + X(i);
    Y_n_1 = Y(i);
end

% Calculate autocorrolation funtion of filtered random process
acf_filtered = autocorr(Y, NumLags=50);
figure(7)
stem(acf_filtered, "filled");

% Calculate power density spectrum of filtered random process
S_y = fft(acf_filtered);
figure(8);
fshift = -25:25;
stem(fshift, abs(fftshift(S_y)), "filled");

%% Task five

% Generate X sequence
Binary_X = ones(10^8,1);

% Randomly select half of them. Since P(x=1) = P(x=-1) = 0.5
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
Error = sum(abs(X_hat - Binary_X));
Error_per = Error/10^8;
disp("Error Percentage is:")
disp(Error_per)






