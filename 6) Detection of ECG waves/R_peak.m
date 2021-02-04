function [R_index] = R_peak(signal,Fs)

    b_lp = [1 zeros(1,5) -2 zeros(1,5) 1]/32;
    a_lp = [1 -2 1];

% High pass filter 
    b_hp = [-1/32 zeros(1,15) 1 -1 zeros(1,14) 1/32];
    a_hp = [1 -1];


% Bandpass
    b_cascade = conv(b_lp,b_hp);
    a_cascade = conv(a_lp,a_hp);

    ecg_bp = filter(b_cascade,a_cascade,signal);

%Derivative based filter
b_d = [2 1 0 -1 -2]/8;
a_d = [1];

%Filtering the signal with the derivative based filter 
ecg_filt = filter(b_d,a_d,ecg_bp);


%Squaring the signal 
ecg_squared = ecg_filt.^2;

%Moving window integrator 
b_int = [ones(1,30)]/30;
a_int = [1];

%Filtering with Moving window 
ecg_int = filter(b_int,a_int,ecg_squared);
ecg_int = ecg_int/max(abs(ecg_int));

%define refractory period to which the algorithm wont search for a peak
ref_p = 0.25;

% Threshold equals to half of maximum outpt of the integrator
threshold = 0.5;

[B] = rmoutliers(ecg_int,'mean');
[ypeak,xpeak] = findpeaks(ecg_int,'MinPeakHeight',threshold*max(B),'MinPeakDistance',ref_p*Fs);
[ypeak,xpeak] = findpeaks(ecg_int,'MinPeakHeight',threshold*mean(ecg_int(xpeak)),'MinPeakDistance',ref_p*Fs);


% Correct for R postion in QRS (Delay = delay of all applied filters (Group
% delay = 0,19s (38 samples)

delay = 38;

%Sample delay may not exactly reflect the actual delay as trapezoid leg
%may not have a sharp transition to its base (hard to say which peak of the integrator is a valid peak) --> need to find the optimum
%delay for the R peaks (delay that will lead to the highest amplitude 

%Inital delay is the peak detected from findpeaks function corrected for
%groupdelay of filters 
initial_delay = xpeak-delay; 

%Create a matrix to find maximal points for all detected peaks
%in the original signal (20 samples after the inital delay which
%corresponds to the length of the QRS complex (100ms)
samples =25;

mat = zeros(length(initial_delay),1);
sample_array = [-samples:samples];

for i=1:length(initial_delay)
    [~,delay_index] = max(signal(initial_delay(i)-samples:initial_delay(i)+samples));
    added_delay = sample_array(delay_index);
    mat(i,:)=added_delay;
end

R_index = initial_delay+mat;
end

