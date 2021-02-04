function [x_lag,corr_res,hz,CSD_T,amplitude] = correlation_T11(x1,x2,fs,lag)

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
 

 
 %Initalize reference vector for correlation
 x1_corr = [zeros(1,maxlag) x1 zeros(1,maxlag)];
 
 %initialize vector to store correlation results 
 
 corr_res = zeros(1,maxlag*2+1);
 
 % correlation loop
 for i=1:length(corr_res)
     matrix= x1_corr(i:i+length(x1)-1).*x2;
     bias = sum(matrix~=0);
     if bias==0
         bias=1;
     end
     corr_res(i) = (sum(matrix))/bias;
 end
 
 % Multiply Correlation result by bias to obtain unbiased estimate 
 %corr_res = corr_res.*bias;
 

 % If  specifies lag is 0, then the output will be defined over the
 % maximum lag 
 
 defined_lag = lag*fs;
 
 if defined_lag ==0 || lag == maxlag
     defined_lag = maxlag;
     x_lag = [-defined_lag:defined_lag];
 end
 
 
 % correct correlation for defined lag 
if maxlag*fs>defined_lag
    mid_index= (length(corr_res)+1)/2;
    corr_res = corr_res(mid_index-defined_lag:mid_index+defined_lag);
    x_lag = [-defined_lag:defined_lag]./fs;
elseif maxlag*fs<defined_lag
    n_zeros = ((defined_lag*2)+1-length(corr_res))/2;
    corr_res = [zeros(1,n_zeros) corr_res zeros(1,n_zeros)];
    x_lag = [-defined_lag:defined_lag]./fs;
else
end

a = corr_res;
CSD_T = fft(a);
% Define x_axis for correlation



slen = length(corr_res);
hz= 0:fs/slen:fs-fs/slen;
amplitude = 2*(abs(fft(corr_res))/slen); %normalization of FFT
amplitude = amplitude(1:length(hz));

end


