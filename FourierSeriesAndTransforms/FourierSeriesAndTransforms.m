%% Fourier Series And Transforms
%% 1. Introduction to DFT and FFT (applicable approach)

% In MATLAB both frequency and time axis are discrete. So we should use an
% operation that is applicable to discrete signals. Just fourier series of
% discrete signals has this property.
imshow("./images/DFT_MATLAB.png");

%% 2. Fourier Series of Discrete Time Signals
N = 1e1;
n = -N:N;
N0 = 5;
M0 = 3;
w0 = 2 * pi / N0;
w = M0 * w0;

x = sin(w * n);

figure('Name', 'Discrete Time Signal');
stem(n, x, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Amplitude');
title('Discrete Time Signal');
grid on;

a_k = fft(x(1:N0)) / N0;
k = n(1:N0);
figure('Name', 'Fourier Series of Discrete Time Signal');
stem(k, imag(a_k), 'LineWidth', 1.5);
xlabel('k');
ylabel('Amplitude');
title('Fourier Series of Discrete Time Signal');
grid on;

% Reconstructing the signal from the Fourier series
kernel = exp(1j * w0 * k' * n);
x_reconstructed = a_k * kernel;

figure('Name', 'Reconstructed Discrete Time Signal');
subplot(2, 1, 1);
stem(n, x, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Amplitude');
title('Original Discrete Time Signal');
grid on;
subplot(2, 1, 2);
stem(n, real(x_reconstructed), 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Amplitude');
title('Reconstructed Discrete Time Signal');
grid on;

%% 3. Fourier Transform of Discrete Time Signals (DTFT)

% The Fourier transform of a discrete time signal is a periodic signal with
% the period of 2 * pi. The Fourier transform of a discrete time signal is
% defined in the interval -pi to pi.

N = 1e3;
n = 0:N;

a = 0.5;
x = a .^ n;

figure('Name', 'Discrete Time Signal');
stem(n, x, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Amplitude');
title('Discrete Time Signal');
grid on;
xlim([-10, 10]);

w_axis = linspace(-pi, pi, 1e3);
DTFT_x = x * exp(-1j * n' * w_axis);

figure('Name', 'DTFT of Discrete Time Signal');
subplot(2, 1, 1);
plot(w_axis, abs(DTFT_x), 'LineWidth', 1.5);
xlabel('Frequency (rad/s)');
ylabel('Amplitude');
title('Absolute Value of DTFT');
xlim([-pi pi]);
xticks([-pi, -pi / 2, 0, pi / 2, pi]);
xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
grid on;

subplot(2, 1, 2);
plot(w_axis, angle(DTFT_x), 'LineWidth', 1.5);
xlabel('Frequency (rad/s)');
ylabel('Phase');
xlim([-pi pi]);
xticks([-pi, -pi / 2, 0, pi / 2, pi]);
xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
title('Phase of DTFT');
grid on;

%% 4. Fourier Transform of Continuous Time Signals
% Considerations:
% 1. Both positive and negative part of time axis should be considered.
% 2. Nyquist theorem should be considered.
% 3. The main part of the signal should be included in the time axis.
% (For periodic signals, the main part is the period of the signal.)

% Increasing the number of samples -> Detecting higher frequencies
% Increasing time axis interval -> higher frequency resolution

% The amplitude of the frequency spectrum is not that important.

% Loweest non-zero frequency -> period = interval of time axis

fs = 1e3;
t = -1:1 / fs:1;
f0 = 10;
x = cos(2 * pi * f0 * t);

FT_x = fft(x);
f_axis = linspace(-fs / 2, fs / 2, length(FT_x));

figure('Name', 'Fourier Transform of Continuous Time Signal');
plot(f_axis, abs(FT_x), 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Fourier Transform of Continuous Time Signal');
grid on;

% FFT formula  in matlab is not the same as the definition of fft formula; It is the same as
% the definition of fftshift. so we should use fftshift to shift the
% frequency axis to the center.
% DFT of a signal is a periodic signal with the period of fs.
% fft function in matlab calculates the fft in the interval 0 to fs.
% fftshift function calculates the fft in the interval -fs / 2 to fs / 2.

% For more accurate results, we should use a longer time axis.
fs = 1e3;
N = 1e3;
t = -N:1 / fs:N;
f0 = 10;
x = cos(2 * pi * f0 * t);

FT_x = (1 / fs) * fftshift(fft(x));
f_axis = linspace(-fs / 2, fs / 2, length(FT_x));

figure('Name', 'Fourier Transform of Continuous Time Signal');
plot(f_axis, abs(FT_x), 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Fourier Transform of Continuous Time Signal');
xlim([-20, 20]);
grid on;

% Question: Why the delta function in frequency have different maximum value?

%% 5. Fourier Series of Continuous Time Signals
fs = 1e4;
T_eq = 1e2;
t = -T_eq:1 / fs:T_eq;
T0 = 0.5;
duty_cycle = 0.25;
square_wave = square(2 * pi * t / T0, duty_cycle * 100);

FT_square_wave = fftshift(fft(square_wave(1:fs * T0))) / (fs * T0);

m0 = 1 + floor((fs * T0) / 2);
N = 1e2;
a_n = FT_square_wave(m0 - N:m0 + N);
n = -N:N;
figure('Name', 'Fourier Series of Square Wave');
stem(n, real(a_n), 'LineWidth', 1.5);
xlabel('n');
ylabel('Amplitude');
title('Fourier Series of Square Wave');
grid on;

% Reconstructing the square wave from the Fourier series
kernel = exp(1j * 2 * pi * n' * t / T0);
x = a_n * kernel;

figure('Name', 'Reconstructed Square Wave');
plot(t, square_wave, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Amplitude');
title('Reconstructed Square Wave');
grid on;
hold on;
plot(t, real(x), 'LineWidth', 1.5);
legend('Original Square Wave', 'Reconstructed Square Wave');
xlim([-0.5, 0.5]);

%% CTFT  extra points
% CTFT of delta function
fs = 1e3;
N = 100;
t = -N:1 / fs:N;
x = zeros(size(t));
x(t == 0) = fs;
FT_x = fftshift(fft(x)) / fs;
f_axis = linspace(-fs / 2, fs / 2, length(FT_x));

figure('Name', 'CTFT of Delta Function  Wrong Answer');
subplot(221);
plot(f_axis, real(FT_x), 'LineWidth', 1.5);
title('Real part of CTFT of Delta Function');
xlabel('Frequency (Hz)');
ylabel('Real(CTFT)');
grid on;

subplot(222);
plot(f_axis, imag(FT_x), 'LineWidth', 1.5);
title('Imaginary part of CTFT of Delta Function');
xlabel('Frequency (Hz)');
ylabel('Imaginary(CTFT)');
grid on;

subplot(223);
plot(f_axis, abs(FT_x), 'LineWidth', 1.5);
title('Magnitude of CTFT of Delta Function');
xlabel('Frequency (Hz)');
ylabel('Magnitude(CTFT)');
grid on;

subplot(224);
plot(f_axis, angle(FT_x), 'LineWidth', 1.5);
title('Phase of CTFT of Delta Function');
xlabel('Frequency (Hz)');
ylabel('Phase(CTFT)');
grid on;

% Constructing delta function using sinc function
fs = 1e3;
N = 100;
t = -N:1 / fs:N;
x = fs * sinc(fs * t);

FT_x = fftshift(fft(x)) / fs;
f_axis = linspace(-fs / 2, fs / 2, length(FT_x));

figure('Name', 'CTFT of Delta Function(kernel = sinc) Wrong Answer');
subplot(221);
plot(f_axis, real(FT_x), 'LineWidth', 1.5);
title('Real part of CTFT of Delta Function(kernel = sinc)');
xlabel('Frequency (Hz)');
ylabel('Real(CTFT)');
grid on;

subplot(222);
plot(f_axis, imag(FT_x), 'LineWidth', 1.5);
title('Imaginary part of CTFT of Delta Function(kernel = sinc)');
xlabel('Frequency (Hz)');
ylabel('Imaginary(CTFT)');
grid on;

subplot(223);
plot(f_axis, abs(FT_x), 'LineWidth', 1.5);
title('Magnitude of CTFT of Delta Function(kernel = sinc)');
xlabel('Frequency (Hz)');
ylabel('Magnitude(CTFT)');
grid on;

subplot(224);
plot(f_axis, angle(FT_x), 'LineWidth', 1.5);
title('Phase of CTFT of Delta Function(kernel = sinc)');
xlabel('Frequency (Hz)');
ylabel('Phase(CTFT)');
grid on;

% True CTFT of delta function
fs = 1e3;
N = 100;
t = -N:1 / fs:N;

x = zeros(size(t));
x(t == 0) = fs;
FT_x = fftshift(fft(circshift(x, -floor(length(t) / 2)))) / fs;

f_axis = linspace(-fs / 2, fs / 2, length(FT_x));

figure('Name', 'True CTFT of Delta Function');
subplot(221);
plot(f_axis, real(FT_x), 'LineWidth', 1.5);
title('Real part of True CTFT of Delta Function');
xlabel('Frequency (Hz)');
ylabel('Real(CTFT)');
grid on;

subplot(222);
plot(f_axis, imag(FT_x), 'LineWidth', 1.5);
title('Imaginary part of True CTFT of Delta Function');
xlabel('Frequency (Hz)');
ylabel('Imaginary(CTFT)');
grid on;

subplot(223);
plot(f_axis, abs(FT_x), 'LineWidth', 1.5);
title('Magnitude of True CTFT of Delta Function');
xlabel('Frequency (Hz)');
ylabel('Magnitude(CTFT)');
grid on;

subplot(224);
plot(f_axis, angle(FT_x), 'LineWidth', 1.5);
title('Phase of True CTFT of Delta Function');
xlabel('Frequency (Hz)');
ylabel('Phase(CTFT)');
grid on;

