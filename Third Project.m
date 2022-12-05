clear;
clc;

% Part a
% generate raised cosine with alpha=0
p = rcosdesign(0,6,4, 'sqrt');
stem(p);
hold on

% generate raised cosine with alpha=0
p = rcosdesign(0.5,6,4, 'sqrt');
stem(p)
xlabel('Samples');
ylabel('Amplitude');
legend('roll-off factor = 0', 'roll-off factor = 0.5')
title('Impluse response with alpha=0.5')

% Convolve with channel response
close;
p = rcosdesign(0,6,4, 'sqrt');
c = zeros(1, 25);
c(13:15) = [0.5, 1/sqrt(2), 0.5];
h = conv(c, p,'same');
stem(h)
  
% Frequency response
close;
H = fft(h);
n = length(h);
fshift = ((-n/2)+0.5:n/2-0.5);
H = fftshift(H);
stem(fshift, abs(H))
xlabel('f')
ylabel('magnitude')
title('Frequency response of h(t)')
 
% Part c
% composite channele impulse response
close;
p = rcosdesign(0.5,6,4, 'sqrt');
c = zeros(1, 25);
c(13:15) = [0.5, 1/sqrt(2), 0.5];
h = conv(c, p,'same');
matched_filter = fliplr(h);
g = conv(h, matched_filter, 'same');
stem(g)
xlabel('Samples');
ylabel('Amplitude');
title('composite channele impulse response')
 
  
% composite channel frequency response
close;
G = fft(g);
n = length(g);
fshift = ((-n/2)+0.5:n/2-0.5);
G = fftshift(G);
stem(fshift, abs(G))
xlabel('f')
ylabel('magnitude')
title('Frequency response of g(t)')
  
% Generate Gaussian noise
close;
N0 = 1;
n = normrnd(0, N0/2, [1, 1000]);
Beta = [0.5, 1/sqrt(2), 0.5];
w = Beta(1)*n + Beta(2)*circshift(n, 1) + Beta(3)*circshift(n, 2);

% Generate signal x[k]
x = randi([0 1], 1, 1000);
x(x==0)=-1;

% Generate r
r = conv(x, Beta, "same") + w;
plot(r)
  
% Detecting the received signal for SNR=0
close;
g0 = g(13);
s = [1, -1];
x_detected = zeros(size(r));

for i=1:numel(r)
    [m, argmin] = min([abs((r(i)/g0)-1)^2, abs((r(i)/g0)+1)^2]);
    x_detected(i) = argmin; 
end

x_detected(x_detected==2)=-1;
[numErrors,ser] = symerr(x,x_detected)
  
% Calculating simulation error for range of SNR:

SNR = 0:2:10;
N0 = (10.^(SNR/10)).^-1;

simulation_error_rate = zeros(1, numel(N0));
for i=1:numel(N0)

    % generate noise
    n = normrnd(0, N0(i)/2, [1, 1000]);
    Beta = [0.5, 1/sqrt(2), 0.5];
    w = Beta(1)*n + Beta(2)*circshift(n, 1) + Beta(3)*circshift(n, 2);
    
    % generate received signal
    r = conv(x, Beta, "same") + w;
    
    x_detected = zeros(size(r));
    for j=1:numel(r)
        [m, argmin] = min([abs((r(j)/g0)-1)^2, abs((r(j)/g0)+1)^2]);
        x_detected(j) = argmin; 
    end
    x_detected(x_detected==2)=-1;
    [numErrors,ser] = symerr(x,x_detected);
    simulation_error_rate(i) = ser;
end

% Calculating theoritical error for range of SNR:
sigma = sqrt(N0./2);
Theoritical_error_rate = qfunc(1./sigma);

% plot the results
close;
semilogy(SNR, simulation_error_rate,'-*')
hold on
semilogy(SNR, Theoritical_error_rate,'-o')
xlabel("SNR");
ylabel("Error rate");
title("Throritical and simulation error rate vs SNR");
legend("simluation error rate", "ideal BPSK error rate");

