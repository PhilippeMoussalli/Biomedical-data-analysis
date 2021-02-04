function [pxx,f,k] = PSDestimation(signal,window,overlap, w_type, fs)
    %   Inputs:
    %        -  signal  : Input signal
    %        -  window  : Window length in seconds
    %        -  overlap : Overlap in proportion (overlap/window)
    %        -  w_type  : Window type as string (select between 'hamming','hanning','bartlett','rectwin').
    %   Outputs:
    %        -  pxx     : PSD estimation
    %        -  f       : frequency vector

signal  =   signal-mean(signal)     ;
limit   =   window*fs               ;
k       =   1                       ;

while limit<length(signal)
    segment =   signal(floor(limit-window*fs+1:limit))         ;
    Win     =   eval([w_type '(length(segment));'])     ;
    %%
    %Complete the code here to compute the periodogram of each segment.
    %Store it in a matrix to then compute an average of the individual
    %estimations:
    seglen = length(segment);
    power = abs(mean(Win.^2));                %Average power of the window
    applied_window =(segment.*Win);
    xdft = fft(applied_window);
    %xdft = xdft(1:floor(seglen/2+1));
    [pxx(:,k)] = (1/seglen*power)*abs(xdft).^2;    %   Compute the periodogram for the k-th segment 
    k       =   k+1                                     ;
    limit   =   limit+round(overlap*window*fs)          ;
end

pxx = mean(pxx,2);
pxx = pxx(1:floor(seglen/2+1)); 
pxx(2:end-1) = 2*pxx(2:end-1);
           
f = 0:fs/seglen:fs/2;                          %   Calculate the vector f