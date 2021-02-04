
clear; close all; clc; 
addpath ./FastICA_25/

%% Initialization

% --- Load the signal in data_ex5_1_2.mat

load('EEG_eyeblink_data.mat');

%% Vizualization of the data

% --- Plot the EEG data using eegPlot
n=length(eeg);
t=linspace(0,n/fs,n);

eegPlot(eeg,channels)

%% Demixing the EEG and identifying the eye blinks

% --- Demix the EEG using fastICA

[Y,A,W,mm]= fastica(eeg);
 Y = zscore(Y,0,2);
 
% --- Plot the temporal information using eegPlot
figure
eegPlot(Y) 

% --- Plot the spatial/topographic information using topoPlot
figure
topoPlot(A)
A([3],1)=0;
% --- Select the eye blink component based on the temporal and spatial
% information

%Given that the eye blink is more prominent in fpz channel, we will select
%that signal to match it with the detected components of ICA to locate
%which component corresponds to the eyeblink

[n_components_detected,~] = size(Y); %ICA does not always output 12 components  

fpz_idx = find(contains(channels,'Fpz'));
denoise_mat = zeros(n_components_detected,1);

for i=1:length(denoise_mat) 
    
      denoise_mat(i,1) = corr(Y(i,:)',eeg(fpz_idx,:)');
   
end

[~,idx_noise] = max(abs(denoise_mat));

%% Removal of the eye blink artifacts

% --- Remove the selected artifact component(s)

A(:,idx_noise)=0;

% --- Reconstruct the EEG without artifacts

eeg_filt= A*Y;

% --- Plot the cleaned reconstructed EEG

figure
eegPlot(eeg_filt,channels) 