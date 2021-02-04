function [amplitude,psd] = FourierT(x,fs)  %Returns amplitude response and PSD in dB
%Fourier, input 
slen = length(x);
hz= linspace(0, fs/2,slen/2);
amplitude = abs(fft(x))/slen; %normalization of FFT
amplitude= amplitude(1:length(hz));
psd= 10*log10((abs(amplitude).^2));
end

%% Fourier plot 
% figure
% plot(hz,amplitude)
% xlim([2 max(hz)])
% grid on
% title('FFT')
% xlabel('Frequency (Hz)')
% ylabel('Amplitude (a.u)')

%% psd plot 
% figure('Name','Power spectrum','NumberTitle','off');
% plot(hz,psd)
% grid on
% title('Power spectrum')
% xlabel('Frequency (Hz)')
% ylabel('Power/Frequency (dB/Hz)')
% axis tight;
     