function [psdx,hz] = psd_W(x,n,fs)

xdft = fft(x,n);
xdft = xdft(1:floor(200/2+1));
psdx = (1/(fs*length(x)) * abs(xdft).^2);
psdx(2:end-1) = 2*psdx(2:end-1);
hz = fs*(0:(n/2))/n;
end

