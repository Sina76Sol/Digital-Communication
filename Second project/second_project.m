% Second Project
% Student name: Reza Soleimani Majd
% Student ID: 100880218
clear;
clc;

%% Task 1
% Performance of QPSK and 4PAM constellations

Symbols = ["00", "01", "10", "11"];

% Creating symbols with equal probabilities:
% ******************* Change 3 to 7 later ********************
X1 = repmat([sqrt(2)/2, sqrt(2)/2], 2.5*10^3, 1);
X2 = repmat([-sqrt(2)/2, sqrt(2)/2], 2.5*10^3, 1);
X3 = repmat([-sqrt(2)/2, -sqrt(2)/2], 2.5*10^3, 1);
X4 = repmat([sqrt(2)/2, -sqrt(2)/2], 2.5*10^3, 1);

X = [X1; X2; X3; X4];

% Shuffle the symbols
shuffled_rows = randperm(size(X, 1));
X = X(shuffled_rows, :);

Y = string(zeros(size(X, 1), 1));
for i = 1:size(X, 1)
    Y(i) = detector(X(i,:));
end
