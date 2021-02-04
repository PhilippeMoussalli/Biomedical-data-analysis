function [amplitude,psdx_n,hz,psdx] = FourierT(x,fs)  %Returns amplitude response and PSD in dB
%psdx_n is the psd up to Nyquist for plotting and psdx is the PSD of the
%full spectrum 

%Fourier, input 
slen = length(x);
hz= 0:fs/slen:fs/2;
amplitude = 2*(abs(fft(x))/slen); %normalization of FFT
amplitude = amplitude(1:length(hz));


%PSD
%PSDX_N
xdft = fft(x);
xdft = xdft(1:floor(slen/2+1)); 
psdx_n = (1/(fs*slen)) * abs(xdft).^2;
psdx_n(2:end-1) = 2*psdx_n(2:end-1);

%PSDX
psdx = (1/(fs*slen)) * abs(fft(x)).^2;  %maybe mutliply by 2 
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
     