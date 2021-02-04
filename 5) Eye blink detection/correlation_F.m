function [x_lag,CCF_F,hz,CSD_F] = correlation_F(x,y,fs,lag)

%  Zero padding to maximum possible lag
 len_dif= length(x)-length(y);
 
 % add zeros if arrays are of different sizes
 if len_dif<0
     x = [x zeros(1,abs(len_dif))];
 elseif len_dif>0
     y = [y zeros(1,len_dif)];
 else
 end 
 
maxlag= (length(x))-1;

% Zero padding to get Number of samples in FFT bigger than len(x)+len(y)-1
x=[zeros(1,maxlag) x];
y=[zeros(1,maxlag) y];

 % Define the maximum lag the correlation function can take 

defined_lag = lag*fs;
if defined_lag ==0 || lag == maxlag
     defined_lag = maxlag;
     x_lag = [-defined_lag:defined_lag];
end
 

n = length(x);
hz = 0:fs/n:fs-fs/n;
amp1 = (fft(x));
amp2= (fft(y));
CSD_F = ((amp1).*conj(amp2));
CCF_F1 = ifft(CSD_F);

% Correct indices location based on defined lag 
CCF_F = [CCF_F1(maxlag+2:end) CCF_F1(1:maxlag+1)];

if maxlag*fs>defined_lag
    mid_index= (length(CCF_F)+1)/2;
    CCF_F = CCF_F(mid_index-defined_lag:mid_index+defined_lag);
     x_lag = [-defined_lag:defined_lag]./fs;
elseif maxlag*fs<defined_lag
    n_zeros = ((defined_lag*2)+1-length(CCF_F))/2;
    CCF_F = [zeros(1,n_zeros) CCF_F zeros(1,n_zeros)];
     x_lag = [-defined_lag:defined_lag]./fs;
else
end

end

%%2.2
% 
% [corr,lag] = xcorr([1,2,3],[1,2,3],'unbiased');
% figure
% plot(lag,corr)
% a= [ zeros(1,3) 1 2 3 5];
% b = [zeros(1,3) 3 2 1 4];
% a1= [1 2 3 5];
% b1 = [3 2 1 4];
% [test,lag]= xcorr(a1,b1);
% af = fft(a);
% bf= fft(b);
% af1 = fft(a1);
% bf1= fft(b1);
% mult = af.*conj(bf);
% mult1 = af1.*conj(bf1);
% test = test
% %nozero = ifft(mult1)
% zero =ifft(mult)