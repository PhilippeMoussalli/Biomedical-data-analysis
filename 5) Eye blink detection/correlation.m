function [x_lag,corr_res,hz,CSD,amplitude,CSD2] = correlation(x1,x2,fs,lag)

 len_dif= length(x1)-length(x2);
 
 % add zeros of arrays are of different sizes
 if len_dif<0
     x1 = [x1 zeros(1,abs(len_dif))];
 elseif len_dif>0
     x2 = [x2 zeros(1,len_dif)];
 else
 end 
 
 % Define the maximum lag the correlation function can take 
 maxlag= (length(x1))-1;
 
 %Initialize N and k for unbiased estimation (Normalization)
 
 N = length(x1);
 k = abs([-maxlag:maxlag]);
 bias= 1./(N-k);
 
 %Initalize reference vector for correlation
 x1_corr = [zeros(1,maxlag) x1 zeros(1,maxlag)];
 
 %initialize vector to store correlation results 
 
 corr_res = zeros(1,maxlag*2+1);
 
 % correlation loop
 for i=1:length(corr_res)
     corr_res(i)= sum(x1_corr(i:i+length(x1)-1).*x2);
 end
 
 % Multiply Correlation result by bias to obtain unbiased estimate 
 corr_res = corr_res.*bias;
 
 
 % If  specifies lag is 0, then the output will be defined over the
 % maximum lag 
 
 defined_lag = lag*fs;
 
 if defined_lag ==0
     defined_lag = maxlag;
 end
 
 
 % correct correlation for defined lag 
if maxlag>defined_lag
    mid_index= (length(corr_res)+1)/2;
    corr_res = corr_res(mid_index-defined_lag:mid_index+defined_lag);
    
elseif maxlag<defined_lag
    n_zeros = ((defined_lag*2)+1-length(corr_res))/2;
    corr_res = [zeros(1,n_zeros) corr_res zeros(1,n_zeros)];
    
else
end

% Define x_axis for correlation
x_lag = [-defined_lag:defined_lag]./fs;


slen = length(corr_res);
hz= 0:fs/slen:fs/2;
amplitude = 2*(abs(fft(corr_res))/slen); %normalization of FFT
amplitude = amplitude(1:length(hz));

%CSD
xdft = fft(corr_res);
xdft = xdft(1:floor(slen/2+1));
CSD = (1/(fs*slen)) * abs(xdft).^2;
CSD(2:end-1) = 2*CSD(2:end-1);
CSD2 = fft(corr_res);
end


