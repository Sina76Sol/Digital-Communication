c = [sqrt(2)/2+sqrt(2)/2i, -sqrt(2)/2+sqrt(2)/2i,...
    -sqrt(2)/2-sqrt(2)/2i, sqrt(2)/2-sqrt(2)/2i];
sigpower = pow2db(mean(abs(c).^2));
M = length(c);

data = randi([0 M-1],10^4,1);
modData = genqammod(data,c);

rxSig = awgn(modData,3,sigpower);

% h = scatterplot(rxSig);
% hold on
% scatterplot(c,[],[],'r*',h)
% grid

%% Experimental error rate
SNR = 0:3:15;
error_rate = zeros(size(SNR));
for i = 1:length(SNR)
    noisy_sig = awgn(modData,SNR(i),sigpower);
    detected = genqamdemod(noisy_sig,c);
    [numErrors,ser] = symerr(data,detected);
    error_rate(i) = ser;
end

plot(error_rate)
%% Imperical error rate
delta =sqrt(2);
Bit_energy = 0.5;
sigma = sqrt(Bit_energy./(10.^(SNR./10)));
Pe = 1-((1-qfunc(delta./(2*sigma)).^2))

