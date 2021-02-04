

%% Clear all previous results and load data

clear all % Clear variables
close all % Close figures
clc

load('ecg2.mat')
resp = val(1,:);
ecg = val(2,:);
fs = 250;

% Example : plot of the original ECG and respiration signal

slen = length(ecg);
t = [1:slen]/fs;
figure
plot(t, ecg)
hold on
plot(t, resp)
xlabel('Time in seconds');
ylabel('Signal (a.u.)');
legend('ECG', 'Respiration')
axis tight;
hold off 
% Plot the power spectrum of the original signal to verify which notch
% filter is needed : let's use fft() to calculate the power
% spectrum of the ECG signal
hz= linspace(0, fs/2,slen/2);

amplitude = abs(fft(ecg))*2; %normalization of FFT
amplitude= amplitude(1:length(hz));
psd_Original = (1/(slen)) * abs(amplitude).^2;
figure
plot(hz,10*log10(psd_Original))
figure('Name','Power spectrum','NumberTitle','off');
plot(hz,10*log10(psd_Original))
grid on
title('Periodogram Using FFT')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
axis tight;
% Power line frequency (notch frequency) was observed at 60 hz 
%% Design the notch filter and apply it to the ECG signal. Observe the effect in the time and frequency domain

Notch_angle = 2*pi*(60/fs);
z1 = cos(Notch_angle)+ i*sin(Notch_angle);
z2 = cos(Notch_angle)- i*sin(Notch_angle);

DC_response = 1+(-z1-z2) +1;

b = [1 -z1-z2  z1*z2]/DC_response;  % Compensation for DC gain 
a = [1];

[h,w]=freqz(b,a,floor(slen/2));
figure('Name','Frequency response of notch filter with zeroes only','NumberTitle','off');

subplot(2,1,1)
plot(hz(1:length(hz)),20*log10(abs(h)));
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('frequency response plot')
axis([0 max(hz) -30 10])

subplot(2,1,2)
plot(hz, 360/(2*pi)*angle(h))
xlabel('Frequency (Hz)')
ylabel('Phase (Degrees)')
title('frequency response plot')
axis([0 max(hz) -120 120])
% Filter ECG with notch filter (only zeros)

y = filter(b,a,ecg);


%% Introduce poles in the notch filter. Try at least 3 values between 0.8 and 0.99 and study the effect. Which filter gives the best result

% Design and plot of notch filters with different radis
a_vector = [0.8 0.95 0.99]

za11= a_vector(1)*z1;
za12= a_vector(1)*z2;

za21= a_vector(2)*z1;
za22= a_vector(2)*z2;

za31= a_vector(3)*z1;
za32= a_vector(3)*z2;

% Compensation for DC gain for different filters (z=1 for H(z)) 
DC_response1 = (1+(-za11-za12)+(za11*za12))/(1+a_vector(1)*(-za11-za12)+ (a_vector(1)^2*(za11*za12)));
DC_response2 = (1+(-za21-za22)+(za21*za22))/(1+a_vector(2)*(-za21-za22)+ (a_vector(2)^2*(za21*za22)));
DC_response3 = (1+(-za31-za32)+(za31*za32))/(1+a_vector(3)*(-za31-za32)+ (a_vector(3)^2*(za31*za32)));


% Initalizing zeros and poles for notch filter with different radis
b_notch = [1 -z1-z2  z1*z2];  


a1 = [1 -za11-za12  za11*za12]*DC_response1;
[h1,w1]=freqz(b_notch,a1,floor(slen/2));
a2 = [1 -za21-za22  za21*za22]*DC_response2;
[h2,w2]=freqz(b_notch,a2,floor(slen/2));
a3 = [1 -za31-za32  za31*za32]*DC_response3;
[h3,w3]=freqz(b_notch,a3,floor(slen/2));

% Filtering with poles
y1 = filter(b_notch,a1,ecg);
y2 = filter(b_notch,a2,ecg);
y3 = filter(b_notch,a3,ecg);


% Frequecny response plotting 
figure('Name','Frequency response of filters','NumberTitle','off');
subplot(2,1,1)
hold on 
plot(hz(1:length(hz)),20*log10(abs(h)));
plot(hz(1:length(hz)),20*log10(abs(h1)));
plot(hz(1:length(hz)),20*log10(abs(h2)));
plot(hz(1:length(hz)),20*log10(abs(h3)));
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('Magnitude')
axis([0 max(hz) -60 10])
legend('zero','Pole at 0.8', 'Pole at 0.95','Pole at 0.99')
hold off 

% Phase response plotting 
subplot(2,1,2)
hold on 
plot(hz, 360/(2*pi)*angle(h))
plot(hz, 360/(2*pi)*angle(h1))
plot(hz, 360/(2*pi)*angle(h2))
plot(hz, 360/(2*pi)*angle(h3))
xlabel('Frequency (Hz)')
ylabel('Phase (Degrees)')
title('Phase')
legend('zero','Pole at 0.8', 'Pole at 0.95','Pole at 0.99')
axis([0 max(hz) -120 120])
hold off 

% Plotting ECG signal before and after filtering  

figure('Name','ECG signal before and after filtering (zeros only and 0.8 radius)','NumberTitle','off');

subplot(2,1,1)
hold on 
plot(t(3:end), ecg(3:end))
axis tight;
plot(t(3:end), y(3:end))
xlabel('Time in seconds');
ylabel('Signal (a.u.)');
title('Original ECG signal filtered with only zeros')
legend('Original signal','zero')
axis tight;
hold off 

subplot(2,1,2)
hold on 
plot(t(3:end), ecg(3:end))
plot(t(3:end), y1(3:end))
xlabel('Time in seconds');
ylabel('Signal (a.u.)');
title('Original ECG signal filtered with pole at 0.8')
legend('Original signal','Pole at 0.8')
axis tight;
hold off 

figure('Name','ECG signal before and after filtering (0.8 and 0.95 radius)','NumberTitle','off');

subplot(2,1,1)
hold on 
plot(t(3:end), ecg(3:end))
plot(t(3:end), y2(3:end))
xlabel('Time in seconds');
ylabel('Signal (a.u.)');
title('Original ECG signal filtered with pole at 0.95')
legend('Original signal','Pole at 0.95')
axis tight;
hold off 

subplot(2,1,2)
hold on 
plot(t(3:end), ecg(3:end))
plot(t(3:end), y3(3:end))
xlabel('Time in seconds');
ylabel('Signal (a.u.)');
title('Original ECG signal filtered with pole at 0.99')
legend('Original signal','Pole at 0.99')
axis tight;
hold off 

% Plotting power spectrum of ECG before and after filtering 

%psd initalization 
amplitude_y = abs(fft(y))*2; 
amplitude_y= amplitude_y(1:length(hz));
amplitude_y1 = abs(fft(y1))*2; 
amplitude_y1= amplitude_y1(1:length(hz));
amplitude_y2 = abs(fft(y2))*2; 
amplitude_y2= amplitude_y2(1:length(hz));
amplitude_y3 = abs(fft(y3))*2; 
amplitude_y3= amplitude_y3(1:length(hz));

psd_y = (1/(slen)) * abs(amplitude_y).^2;
psd_y1 = (1/(slen)) * abs(amplitude_y1).^2;
psd_y2 = (1/(slen)) * abs(amplitude_y2).^2;
psd_y3 = (1/(slen)) * abs(amplitude_y3).^2;


figure('Name','Power spectrum (zeros and 0.8 notch)','NumberTitle','off');
subplot(2,1,1)
hold on 
plot(hz,10*log10(psd_Original))
plot(hz,10*log10(psd_y))
grid on
title('Power spectrum of original signal with filered signal (zeros only)')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
legend('Original signal','Zeros')
axis tight;

subplot(2,1,2)
hold on 
plot(hz,10*log10(psd_Original))
plot(hz,10*log10(psd_y1))
grid on
title('Power spectrum of original signal with filered signal (pole at 0.8)')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
legend('Original signal','Pole at 0.8')
axis tight;


figure('Name','Power spectrum (0.95 and 0.99 notch)','NumberTitle','off');
subplot(2,1,1)
hold on 
plot(hz,10*log10(psd_Original))
plot(hz,10*log10(psd_y2))
grid on
title('Power spectrum of original signal with filered signal (pole at 0.95)')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
legend('Original signal','Pole at 0.95')
axis tight;

subplot(2,1,2)
hold on 
plot(hz,10*log10(psd_Original))
plot(hz,10*log10(psd_y3))
grid on
title('Power spectrum of original signal with filered signal (pole at 0.99)')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
legend('Original signal','Pole at 0.99')
axis tight;

%% Apply the filter to the provided respiration signal.

%PSD of respiration
amplitude_resp = abs(fft(resp))*2; 
amplitude_resp= amplitude_resp(1:length(hz));
psd_resp = (1/(slen)) * abs(amplitude_resp).^2;


figure('Name','PSD resp','NumberTitle','off');
plot(hz,10*log10(psd_resp))
grid on
title('Power spectrum of respiration signal')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
axis tight;

% Filter respiration signal with notch at 0.95 
resp_filtered = filter(b_notch,a2,resp);

% Plots of ECG and respiration before and after filtering 

figure
plot(t, ecg)
hold on
plot(t, y2)
plot(t, resp)
plot(t, resp_filtered)
xlabel('Time in seconds');
ylabel('Signal (a.u.)');
legend('ECG', 'ECG filtered at 0.95','respiration','respiration filtered')
axis tight;
hold off 

%% show overlap betweem PSD of EEG and PSD of repiration 


amplitude_resp1 = abs(fft(resp_filtered))*2; 
amplitude_resp1= amplitude_resp1(1:length(hz));
psd_resp1 = (1/(slen)) * abs(amplitude_resp1).^2;


figure('Name','PSD resp','NumberTitle','off');
hold on 
plot(hz,10*log10(psd_y2))
plot(hz,10*log10(psd_resp1))
grid on
title('Power spectrum of respiration signal')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
axis tight;

