function [amplitude,hz] = FourierT(x,fs)
%Fourier, input 
slen = length(x);
hz= linspace(0, fs/2,slen/2);
amplitude = abs(fft(x))/slen; %normalization of FFT
amplitude= amplitude(1:length(hz));
end

% %
% figure
% plot(hz,amplitude)
% xlim([2 max(hz)])
% grid on
% title('FFT')
% xlabel('Frequency (Hz)')
% ylabel('Amplitude (a.u)')
     