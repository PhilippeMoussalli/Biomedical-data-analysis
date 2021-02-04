---------------------------------------

%% 3.0 Initialization
clearvars;
close('all');

% --- load the signal in _data_ex2_3.mat_

load('Eye blink detection_data.mat')

fs = EEG.fs;
data = EEG.x;
labels = EEG.label;

% --- plot the multi-channel data
% --- (use the function _eegplot2019_)

figure('name','EEG data')
eegplot2019(data,fs,labels)
title('EEG')

% --- extract an eyeblink template from the data
chan_name1 = 'F4'; % pick a certain channel (don't select another template for this exercise)
chan_idx1 = find(cellfun(@(x)~isempty(x),strfind(lower(EEG.label),lower(chan_name1))),1);

blinkstart = 0.9;
template = EEG.x(chan_idx1,round(blinkstart*EEG.fs+(1:0.5*EEG.fs)));
t_template = (1:length(template))/EEG.fs;
lag = (length(t_template)/2)-1; % (Measured in samples) Need to correct for lag as dot product in CCF begins with signal alligned in the middle 

% --- extract the signal at the selected channel

sig1 = EEG.x(chan_idx1,:);
t_axis = (1:length(sig1))/EEG.fs; % construct a time axis for the signal

%% 3.1 Template matching

% --- select a channel from the parietal or occipital lobe (any channel name 
% --- that is a combination of 'P' and/or 'O', respectively) and compute the
% --- normalized/rescaled CCF between the template and 
% --- a) the signal from that channel
% --- b) the signal from the template's channel

chan_name2 = 'O2'; % pick a certain channel (don't select another template for this exercise)
chan_idx2 = find(cellfun(@(x)~isempty(x),strfind(lower(EEG.label),lower(chan_name2))),1);

sig2 = EEG.x(chan_idx2,:);% ... (signal from a particular channel)


[x_lag1,corr_res1,~,~,~] = correlation_T(sig1,template,fs,0);
[x_lag2,corr_res2,~,~,~] = correlation_T(sig2,template,fs,0);

zero_index = find(x_lag1 ==0);
ccf1 = corr_res1(zero_index:end);% ...(the CCF between the template and sig1)
ccf2 = corr_res2(zero_index:end);% ... (the CCF between the template and sig1)

% --- plot the signals and the normalized CCF 
figure
subplot(2,1,1)
hold on 
plot(t_axis,sig1/max(abs(sig1)),'b')
plot(t_axis,ccf1/max(abs(ccf1)),'r')
legend('Signal','Cross correlation')
xlabel('Time')
ylabel('Amplitude (uV)')
title('CCF with sig1 (Channel F4)')
%ylim([-200 200])

subplot(2,1,2)
hold on 
plot(t_axis,sig2/max(abs(sig2)),'b')
plot(t_axis,ccf2/max(abs(ccf2)),'r')
legend('Signal','Cross correlation')
xlabel('Time')
ylabel('Amplitude (uV)')
title('CCF with sig1 (Channel O2)')

%% 2.2 Matched filtering

% --- compute the impulse response h_mf from the matched filter

h_mf = ifft( conj( fft(template) ) );

%  by taking the complex conjugate of the Fourier coefficients,
% the template gets flipped in the time domain

% --- plot h_mf together with the template

figure,hold('on')
plot((1:length(template))/EEG.fs,template)
plot((1:length(template))/EEG.fs,h_mf)
xlabel('Time')
ylabel('Amplitude (uV)')
title('Template and matched filter impulse response')
legend('template','matched filter impulse response','Reversed template')

% --- apply the matched filter to a selected signal
% --- note: h_mf is a FIR filter, hence you can use the _fftfilt_ command

%help('fftfilt')
sig1_mf = fftfilt(h_mf,sig1);

% --- plot the matched filter's output, together with the signal and the CCF

figure,hold('on')
plot(t_axis,sig1/max(abs(sig1)))
plot(t_axis,sig1_mf/max(abs(sig1_mf)))
plot(t_axis,ccf1/max(abs(ccf1)))
% <plot also the normalized CCF between template and sig0>
legend('signal','matched filter output','normalized CCF')
title('Matched filter vs. CCF of a signal')
xlabel('Time')
ylabel('Amplitude (Normalized)')
axis('tight')


compute the onsets of
% the eyeblinks via a thresholding procedure

threshold = 0.5;
[ypeak,xpeak] = findpeaks(sig1_mf,'MinPeakHeight',threshold*max(sig1_mf),'MinPeakDistance',length(template)); %Avoid detecting multile peaks that are close to each other 

% based on the locations of the peaks in the output, find now where the
% _onsets_ (and not the peaks!) of the eyeblinks are

xonset = xpeak-length(template);

% --- superimpose markers at the onsets of the eyeblinks
figure,hold('on')
plot(t_axis,sig1/max(abs(sig1)))
plot(t_axis,sig1_mf/max(abs(sig1_mf)))
stem(t_axis(xonset),0.5*ones(size(xonset)),'-ko','filled')
legend('signal','matched filter output','eyeblink onsets')
xlabel('Time')
ylabel('Amplitude (Normalized)')
title('Template and matched filter impulse response')
title('Eyeblink detection using matched filter')
axis('tight')
