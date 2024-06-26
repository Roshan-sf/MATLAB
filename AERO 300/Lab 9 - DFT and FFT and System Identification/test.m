% Load variables from the .mat file
load('fourierData.mat');

% MATLAB code for DFT and IDFT processing
clc;

% Use the loaded variables: Fs, L, S, T, X, t

% Sampling frequency and number of points
N = L; % Number of points from the loaded variable L

% Perform DFT on the signal X
Xk = calcdft(X, N);

disp('DFT X(k): ');
disp(Xk);

% Calculate magnitude and phase of the DFT
mgXk = abs(Xk);
phaseXk = angle(Xk);

% Frequency vector
k = 0:N-1;

% Plot magnitude and phase of the DFT
figure;
subplot(2,1,1);
stem(k, mgXk);
title('DFT Magnitude');
xlabel('Frequency');
ylabel('Magnitude');

subplot(2,1,2);
stem(k, phaseXk);
title('DFT Phase');
xlabel('Frequency');
ylabel('Phase');

% Perform IDFT on the DFT result
xn = calcidft(Xk);

disp('Inverse DFT x(n): ');
disp(xn);

% Time vector
n = 0:N-1;

% Plot the IDFT result
figure;
stem(n, xn);
title('Inverse DFT');
xlabel('Time');
ylabel('Amplitude');

% Function to calculate DFT
function [Xk] = calcdft(xn, N)
    L = length(xn);
    if N < L
        error('N must be greater than or equal to L!!');
    end
    x1 = [xn, zeros(1, N - L)];
    W = zeros(N, N);
    for k = 0:N-1
        for n = 0:N-1
            W(k+1, n+1) = exp(-1i * 2 * pi * n * k / N);
        end
    end
    disp('Transformation matrix for DFT');
    disp(W);
    Xk = W * x1.'; % Ensure x1 is a column vector for multiplication
end

% Function to calculate IDFT
function [xn] = calcidft(Xk)
    N = length(Xk);
    IT = zeros(N, N);
    for k = 0:N-1
        for n = 0:N-1
            IT(k+1, n+1) = exp(1i * 2 * pi * n * k / N);
        end
    end
    disp('Transformation Matrix for IDFT');
    disp(IT);
    xn = (IT * Xk) / N; % Ensure Xk is a column vector for multiplication
end
