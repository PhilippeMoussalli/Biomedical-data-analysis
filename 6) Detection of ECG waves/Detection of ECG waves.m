

close all
clear 
clc
% Initializing parameters

load('ECG_data.mat')                 ;   % llowing the ECGsignal
%ecg     =   signal;
Fs      = 200                   ;   % Sampling frequency
time    = (0:length(ecg)-1)/Fs  ;   % Creating a Time Vector
n=length(ecg);

figure(1)
plot(time,ecg);
axis 'tight'
grid on
xlabel('time(seconds)')
ylabel('Amplitude [A.U.]')
title('ECG signal for exercise 1')


%% Pan-Tompkins algorithm


% Low-pass Filter

b_lp = [1 zeros(1,5) -2 zeros(1,5) 1]/32;
a_lp = [1 -2 1];

% High pass filter 
b_hp = [-1/32 zeros(1,15) 1 -1 zeros(1,14) 1/32];
a_hp = [1 -1];


% Bandpass
b_cascade = conv(b_lp,b_hp);
a_cascade = conv(a_lp,a_hp);

ecg_bp = filter(b_cascade,a_cascade,ecg);

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
[ypeak,xpeak] = findpeaks(ecg_int,'MinPeakHeight',threshold*max(ecg_int),'MinPeakDistance',ref_p*Fs);

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
samples =20;

mat = zeros(length(initial_delay),1);

for i=1:length(initial_delay)
    [~,delay_index] = max(ecg(initial_delay(i):initial_delay(i)+samples-1));
    mat(i,:)=delay_index;
end

R_index = initial_delay+mat-1;

%% Find Q 

matQ = zeros(length(R_index),1);

% Find minimum value 60 ms after R peak which correspond to Q 

for i=1:length(R_index)
    [~,Q_index] = min(ecg(R_index(i)-(0.06*Fs):R_index(i)));
    matQ(i,:)=Q_index;
end

% Define Q index 
interval_length = length(ecg(R_index(1)-(0.06*Fs):R_index(1))); %Correct for reversed interval 
Q_index = (R_index-interval_length)+ matQ; %Correction to get proper index 


%% Find S 

% Use the output of the derivative filter to localize 'S' (Corresponds to
% the lowest peak after R)

% Easier to find it from the derivative response 
[ymin_S,xmin_S] = findpeaks(-ecg_filt,'MinPeakHeight',threshold*max(-ecg_filt),'MinPeakDistance',ref_p*Fs); %Avoid detecting multile peaks that are close to each other 

%Correct for delay of 23 samples (115 msec) (Bandpass+Differentiator) 
delay_S = 23;
inital_S = xmin_S-delay_S;

% Optimize to find minimum value in an interval of 5 samples
samples_S =5;

mat_S = zeros(length(inital_S),1);

for i=1:length(inital_S)
    [~,delay_index_S] = min(ecg(inital_S(i):inital_S(i)+samples_S-1));
    mat_S(i,:)=delay_index_S;
end

S_index = inital_S+mat_S-1;

%% P and T wave detection 

% 1) Delete QRS complex by taking 5 ms before R peak and 7 ms after. The
% baseline as the average signal taken between 10 ms and 5 ms and  before the R peak.

%% T wave detection 

ecg2=ecg; %Create a copy of ECG to manipulare 
ecg2=ecg2+abs(ecg(1)); %Shift ECG baseline to 0

% Replace QRS with 0 
for i=1:length(R_index)
    ecg2(Q_index(i)-5:S_index(i)+5) =0;  %+5 -5 to account for the full range of Q and S
end

%Pass signal through derivative based filter,square and calculate integral
%(Similar to Gritazli et.al method)
ecg2_filt = filter(b_d,a_d,ecg2);
ecg2_filt_s = ecg2_filt.^2;
ecg_int2 = filter(b_int,a_int,ecg2_filt_s);
ecg_int2 =ecg_int2(1:end-15); %account for undetected R peak
ecg_int2 = ecg_int2/max(abs(ecg_int2));

% Find maximal peak that corresponds to T wave 
[ymax_T,xmax_T] = findpeaks(ecg_int2 ,'MinPeakHeight',threshold*max(ecg_int2),'MinPeakDistance',ref_p*Fs); 

% correct for the delay of the three applied filters 17 samples (85ms)
xmax_T = xmax_T-17;
timeint2 = time(1:length(ecg_int2));  %Time xaxis for Integral result 

%Optimize to obtain the highest T-peak value
samples_T =5;

mat_T = zeros(length(xmax_T),1);

for i=1:length(xmax_T)
    [~,delay_index_T] = max(ecg(xmax_T(i):xmax_T(i)+samples_T-1));
    mat_T(i,:)=delay_index_T;
end

T_index = xmax_T+mat_T-1; %Optimal T index 

% Delete T and replace baseline --> Filter the signal just as in the steps
% used for T detection and then the output of the integral will correspond
% to the peak of the P wave 

%% P wave detection 

%Shift baseline to start around zero 

ecg3=ecg+abs(ecg(1)); 

% Delete QRS and P from the signal 
for i=1:length(R_index)
    ecg3(Q_index(i)-5:T_index(i)+(0.085*Fs)) = 0; 
end

% Specify  window length of the integral filter to have the same width as the P wave (16
% samples~= 80ms)

b_int2 = [ones(1,16)]/16;
a_int2 = [1];

%Calculate delay of integral filter with 16 samples 
[gd]= grpdelay(b_int2,a_int2);
int_delay_16 = gd(1);

% delay for int_filter with 16 samples is 7.5 samples 

%Pass signal through derivative based filter,square and calculate integral 
ecg3_filt = filter(b_d,a_d,ecg3);
ecg3_filt_s = ecg3_filt.^2;
ecg3_filt_sq=sqrt(ecg3_filt_s);
ecg_int3 = filter(b_int2,a_int2,ecg3_filt_s);
ecg_int3 =ecg_int3(1:end-15); %account for undetected R peak
ecg_int3 = ecg_int3/max(abs(ecg_int3));

[ymax_P,xmax_P] = findpeaks(ecg_int3,'MinPeakHeight',threshold*max(ecg_int3),'MinPeakDistance',ref_p*Fs);

%Total delay = 7.5+ 2(differentiator)~= 9 samples
% correct for the delay of the three applied filters (45ms)

P_index = xmax_P-9;


%% Tachogram computation 

R_time=time(R_index);
tachogram = (R_time(2:end)-R_time(1:end-1)); 
%% Figure 1
% A figure showing the steps to detect the R-peaks

figure(2)
ha(1)= subplot(511);
plot(time,ecg);
axis 'tight'
grid on
ylabel('ECG')
xlim([0, 2.2])
title('R peak detection')

ha(2)=subplot(512);
plot(time,ecg_bp);
axis 'tight'
grid on
ylabel('BFP')
xlim([0, 2.2])

ha(3)=subplot(513);
plot(time,ecg_filt);
axis 'tight'
grid on
ylabel('Der')
xlim([0, 2.2])


ha(4)=subplot(514);
plot(time,ecg_squared);
axis 'tight'
grid on
ylabel('Sqr')
xlim([0, 2.2])

ha(5)=subplot(515);
plot(time,ecg_int);
axis 'tight'
grid on
xlabel('Time in seconds')
ylabel('Int')
xlim([0, 2.2])

linkaxes(ha, 'x');
%suptitle('Steps of the Pan-Tomkins algorithm')


% A figure showing the steps to detect the T-peaks

figure(3)
ha(1)= subplot(411);
plot(time,ecg2);
axis 'tight'
grid on
ylabel('ECG')
xlim([0, 2.2])
title('T peak detection')


ha(3)=subplot(412);
plot(time,ecg2_filt);
axis 'tight'
grid on
ylabel('Der')
xlim([0, 2.2])


ha(4)=subplot(413);
plot(time,ecg2_filt_s);
axis 'tight'
grid on
ylabel('Sqr')
xlim([0, 2.2])

ha(5)=subplot(414);
plot(timeint2,ecg_int2);
axis 'tight'
grid on
xlabel('Time in seconds')
ylabel('Int')
xlim([0, 2.2])

linkaxes(ha, 'x');

% A figure showing the steps to detect the P-peaks

figure(4)
ha(1)= subplot(411);
plot(time,ecg3);
axis 'tight'
grid on
ylabel('ECG')
xlim([0, 2.2])
title('P peak detection')


ha(3)=subplot(412);
plot(time,ecg3_filt);
axis 'tight'
grid on
ylabel('Der')
xlim([0, 2.2])


ha(4)=subplot(413);
plot(time,ecg3_filt_s);
axis 'tight'
grid on
ylabel('Sqr')
xlim([0, 2.2])

ha(5)=subplot(414);
plot(timeint2,ecg_int3);
axis 'tight'
grid on
xlabel('Time in seconds')
ylabel('Int')
xlim([0, 2.2])

linkaxes(ha, 'x');
%% Figure 2
% A figure of the ECG signals with the detected R-peaks and the
% corresponding tachogram for the signal

figure(5)
subplot 211
title('ECG signal and its events')
hold on 
plot(time,ecg)  
plot(time(T_index),ecg(T_index),'yo','MarkerFaceColor','y','MarkerEdgeColor','k')
plot(time(P_index),ecg(P_index),'bo','MarkerFaceColor','b','MarkerEdgeColor','k')
plot(time(R_index),ecg(R_index),'ro','MarkerFaceColor','r','MarkerEdgeColor','k')
plot(time(Q_index),ecg(Q_index),'ro','MarkerFaceColor','r','MarkerEdgeColor','k')
plot(time(S_index),ecg(S_index),'ro','MarkerFaceColor','r','MarkerEdgeColor','k')
axis 'tight'
grid on
xlabel('time(seconds)')
ylabel('Amplitude [A.U.]')

legend('ECG signal','T-wave','P-wave','QRS complex')

subplot 212
plot((1:length(tachogram)),tachogram)
axis 'tight'
grid on
xlabel('Time (in second)')
ylabel('Duration (in second)')
title('Tachogram')
%suptitle('Computation of the tachogram in the test signal')