%% EX3.3. Cascaded filters for ECG 

clear all % Clear variables
close all % Close figures
clc
load('ecg3.mat') % Load signals
fs = 1000; % Sampling frequency
signal = ECG23(:,1);
N = length(signal);
hz= linspace(0, fs/2 ,N/2);
t=linspace(1/fs, N/fs, N);

%% Hanning Filter

% Iniitalizing Hanning filter
b_H = [0.25 0.5 0.25];
a_H = [1];

% Mangitude and phase response

[h_H,w_H]=freqz(b_H,a_H,floor(N/2));

% Frequecny response plotting 
figure('Name','Frequency response of Hanning','NumberTitle','off');

subplot(2,1,1)
hold on 
plot(hz(1:length(hz)),20*log10(abs(h_H)));
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('Magnitude')
axis tight
hold off 

% Phase response plotting 
subplot(2,1,2)
hold on 
plot(hz, 360/(2*pi)*angle(h_H))
xlabel('Frequency (Hz)')
ylabel('Phase (Degrees)')
title('Phase')
axis tight
hold off 

% Signal smoothing (Hanning filter)

y_H = filter(b_H,a_H,signal);


% Plot of signals before and after filtering(Hanning filter)
figure('Name','Hanning filter','NumberTitle','off');

subplot(2,1,1)
plot(t, signal)
xlabel('Time in seconds');
ylabel('Signal (a.u.)');
title('Original signal')
axis tight;
subplot(2,1,2)
plot(t, y_H)
xlabel('Time in seconds');
ylabel('Signal (a.u.)');
title('Hanning filter')
axis tight;
hold off 

% Plot of Fourier spectra of signals before and after filtering (Hanning filter)
[A_Signal,psd_signal] = FourierT(signal,fs);
[A_H] = FourierT(y_H,fs);
figure('Name','Hanning filter (Fourier Spectrum)','NumberTitle','off');
hold on 
plot(hz,A_Signal)
plot(hz,A_H)
xlim([2 max(hz)])
grid on
title('FFT')
xlabel('Frequency (Hz)')
ylabel('Amplitude (a.u)')
legend('Original','Hanning')
axis([0 500 0 1])

%% Derivative based filter 

% Iniitalizing of derivative based filter 
b_D = [1 -1];
a_D = [1 -0.99]; % The cutoff frequency is lower for 0.99

% Mangitude and phase response

[h_D ,w_D]=freqz(b_D,a_D,floor(N/2));

% Frequecny response plotting 
figure('Name',' Frequency response of Derivative based filtering','NumberTitle','off');

subplot(2,1,1)
hold on 
plot(hz(1:length(hz)),20*log10(abs(h_D)));
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('Magnitude')
axis tight
hold off 

% Phase response plotting 
subplot(2,1,2)
hold on 
plot(hz, 360/(2*pi)*angle(h_D))
xlabel('Frequency (Hz)')
ylabel('Phase (Degrees)')
title('Phase')
axis tight
hold off 

% Remove baseline drift (Derivative based filter)

y_D = filter(b_D,a_D,signal);

figure('Name','Deriative based filter','NumberTitle','off');

subplot(2,1,1)
plot(t, signal)
xlabel('Time in seconds');
ylabel('Signal (a.u.)');
title('Original signal')
axis tight;
subplot(2,1,2)
plot(t, y_D)
xlabel('Time in seconds');
ylabel('Signal (a.u.)');
title('Derivative based filter')
axis tight;
hold off 

% Plot of Fourier spectra of signals before and after filtering (Derivative based filter)
[A_D] = FourierT(y_D,fs);
figure('Name','Derivative based filtering (Frequency spectrum)','NumberTitle','off');
hold on 
plot(hz,A_Signal)
plot(hz,A_D)
xlim([0 max(hz)])
grid on
title('FFT')
xlabel('Frequency (Hz)')
ylabel('Amplitude (a.u)')
legend('Original','Derivative based filter')
axis([0 5 0 20])

%% Comb filter

% Power spectrum analysis to identify Power line frequency and its harmonics
amplitude = abs(fft(signal))*2; %normalization of FFT
amplitude= amplitude(1:length(hz));
psd= (1/(N)) * abs(amplitude).^2;
figure('Name','Original Signal Power spectrum','NumberTitle','off');
plot(hz,10*log10(psd))
grid on
title('Periodogram Using FFT')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
axis tight;

% Frequencies identifed at (50,150,250,350,450 Hz)

freq_comb = 50; 
comb_mat = zeros(5,2); %Define matrix for transfer function 
coeff_mat = ones(5,3); 

for i_comb = 1:length(comb_mat)
    comb_angle = 2*pi*(freq_comb/1000); 
    comb_mat (i_comb, 1) = cos(comb_angle)+ i*sin(comb_angle);
    comb_mat (i_comb, 2) = cos(comb_angle)- i*sin(comb_angle);
    coeff_mat (i_comb, 2) = -(comb_mat(i_comb, 1))-comb_mat(i_comb, 2);
    freq_comb = freq_comb + 100;
end

%Convolution to obtain polynomial multiplication 
comb_output = (conv(coeff_mat(5,:), conv(coeff_mat(4,:), conv(coeff_mat(3,:), conv(coeff_mat(2,:), coeff_mat(1,:))))));

G = sum(comb_output); %factor to set the gain at DC = 1

b_C = comb_output / G; 
a_C = [1];

% Mangitude and phase response

[h_C,w_C]=freqz(b_C,a_C,floor(N/2));

% Frequecny response plotting 
figure('Name','Frequency response of Comb Filter','NumberTitle','off');

subplot(2,1,1)
hold on 
plot(hz(1:length(hz)),20*log10(abs(h_C)));
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('Magnitude')
axis tight
hold off 

% Phase response plotting 
subplot(2,1,2)
hold on 
plot(hz, 360/(2*pi)*angle(h_C))
xlabel('Frequency (Hz)')
ylabel('Phase (Degrees)')
title('Phase')
axis tight
hold off 

% Eliminating power line and hamronics(Comb filter)

y_C = filter(b_C,a_C,signal);

figure('Name','Comb Filter','NumberTitle','off');

subplot(2,1,1)
plot(t, signal)
xlabel('Time in seconds');
ylabel('Signal (a.u.)');
title('Original signal')
axis tight; 
subplot(2,1,2)
plot(t, y_C)
xlabel('Time in seconds');
ylabel('Signal (a.u.)');
title('Comb filter')
axis tight;
hold off 


% Plot of Fourier spectra of signals before and after filtering (Comb filter)
[A_C] = FourierT(y_C,fs);

figure('Name','Comb Filter (Fourier Spectrum)','NumberTitle','off');
hold on 
plot(hz,A_Signal)
plot(hz,A_C)
xlim([0 max(hz)])
grid on
title('FFT')
xlabel('Frequency (Hz)')
ylabel('Amplitude (a.u)')
legend('Original','Comb filter')
axis([2 max(hz) 0 3])


%% Cascading all the three filters
b_cascade = conv(b_C, conv(b_D, b_H));
a_cascade = conv(a_C, conv(a_D, a_H));

[h_cascade,w_cascade]=freqz(b_cascade,a_cascade,floor(N/2));

% Frequency response plotting of the 3 filters
figure('Name','Frequency response of Cascaded Filters','NumberTitle','off');

subplot(2,1,1)
hold on 
plot(hz(1:length(hz)),20*log10(abs(h_cascade)));
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('Magnitude')
axis tight
hold off 

% Phase response plotting 
subplot(2,1,2)
hold on 
plot(hz, 360/(2*pi)*angle(h_cascade))
xlabel('Frequency (Hz)')
ylabel('Phase (Degrees)')
title('Phase')
axis tight
hold off 

%Plot of cascaded filter
y_cascade = filter (b_cascade, a_cascade, signal);

figure('Name','Cascaded Filter','NumberTitle','off');
hold on
plot(t, signal);
plot(t, y_cascade);
xlabel('Time in seconds');
ylabel('Signal (a.u.)');
title('Use of Cascaded Filter')
axis tight;
legend('Original Signal', 'Filtered Signal');
hold off

%%Overlaid power sepctrum of signal before and after filtering 
[A_Signal_cascaded,psd_signal_cascaded] = FourierT(y_cascade,fs);
figure('Name','PSD spectrum','NumberTitle','off');
hold on 
plot(hz,psd_signal)
plot(hz,psd_signal_cascaded)
grid on
title('Power spectrum')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
legend('Original','Filtered (Cascade)')
axis tight;


%%Group Delay calculation 
figure('Name','Group Delay','NumberTitle','off');
figure
grpdelay(b_cascade, a_cascade);
title('Group Delay Plot')
