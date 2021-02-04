
%% Clear all previous results and load data

clear all % Clear variables
close all % Close figures
clc
load('detectiontask.mat') % Load signals
fs = 250; % Sampling frequency

% Example: Plots the full signal of electrode C3

C3 = EEGdata(5,:); %The signals of individual electrodes are stored as rows in the EEGdata matrix
tt=linspace(1/fs, length(C3)/fs, length(C3)); % Construct time signal
plot(tt,C3)
axis('tight');
xlabel('Time in seconds');
ylabel('Amplitude in \muV');
title('EEG channel - C3');


%% Segment the data for the UL and UR events in channels PO7 and PO8. 
%The segments (or epochs) should start 100 ms before the stimulus and 600 ms after the stimulus. 
%For each epoch, estimate and subtract the baseline, calculated using the mean of the 100ms PRE-stimulus interval.
PO7_index = find(contains(labels,'PO7'));
PO8_index = find(contains(labels,'PO8'));
PO7 = EEGdata(PO7_index,:);
PO8 = EEGdata(PO8_index,:);

%convert UL and UR time instance to time (in ms)

S1_time = (S1*(1/fs))*1000;
S2_time = (S2*(1/fs))*1000;

% create matrix of start and end point of each epoch for UL and UR

S1_time_matrix = zeros(2,length(S1_time));
S2_time_matrix = zeros(2,length(S2_time));

% assign start and end point of epoch for UL and UR

for i=1:length(S1_time)
    S1_time_matrix(1,i)= S1_time(i)-100;
    S1_time_matrix(2,i)= S1_time(i)+600;
end

for i=1:length(S2_time)
    S2_time_matrix(1,i)= S2_time(i)-100;
    S2_time_matrix(2,i)= S2_time(i)+600;
end
    

% Create matrix for each epoch with 3 indices: (Start point,stimulus and
% end point) for S1 and S2

UL_epoch_matrix = zeros(3,length(S1_time));
UR_epoch_matrix = zeros(3,length(S2_time));

for i=1:length(S1_time)
    UL_epoch_matrix(1,i) = (S1_time_matrix(1,i)*fs)/1000;  % index of the start of the pre-stimulus
    UL_epoch_matrix(2,i) = S1(i);                           % index of stimulus
    UL_epoch_matrix(3,i) = (S1_time_matrix(2,i)*fs)/1000;  % index of the end of the pre-stimulus
end


for i=1:length(S2_time)
    UR_epoch_matrix(1,i) = (S2_time_matrix(1,i)*fs)/1000;  % index of the start of the pre-stimulus
    UR_epoch_matrix(2,i) = S2(i);
    UR_epoch_matrix(3,i) = (S2_time_matrix(2,i)*fs)/1000; 
end

% Baseline defintion

UL_baseline = zeros(2,length(S1_time));  % 2 rows (1st for Channel PO7  second for Channel PO8)
UR_baseline = zeros(2,length(S2_time));

% S1 for PO7 and PO8
for i=1:length(S1_time)
    UL_baseline(1,i) =  mean(PO7(UL_epoch_matrix(1,i):UL_epoch_matrix(2,i)));
    UL_baseline(2,i) =  mean(PO8(UL_epoch_matrix(1,i):UL_epoch_matrix(2,i)));
end

% S2 for PO7 and PO8

for i=1:length(S2_time)
    UR_baseline(1,i) =  mean(PO7(UR_epoch_matrix(1,i):UR_epoch_matrix(2,i)));
    UR_baseline(2,i) =  mean(PO8(UR_epoch_matrix(1,i):UR_epoch_matrix(2,i)));
end

% Substact basline and define different epochs for UL and UR 
epoch_index_size = length(PO7(UL_epoch_matrix(1,1):UL_epoch_matrix(3,1)));
PO7_UL_epoch_values = zeros(epoch_index_size,length(S1_time));
PO7_UR_epoch_values = zeros(epoch_index_size,length(S2_time));
PO8_UL_epoch_values = zeros(epoch_index_size,length(S1_time));
PO8_UR_epoch_values = zeros(epoch_index_size,length(S2_time));

PO7_UL_epoch_values1 = zeros(epoch_index_size,length(S1_time));
PO7_UR_epoch_values1 = zeros(epoch_index_size,length(S2_time));
PO8_UL_epoch_values1 = zeros(epoch_index_size,length(S1_time));
PO8_UR_epoch_values1 = zeros(epoch_index_size,length(S2_time));

for i=1:length(S1_time)
    PO7_UL_epoch_values(:,i) = PO7(floor(UL_epoch_matrix(1,i)):floor(UL_epoch_matrix(3,i)))-UL_baseline(1,i);
end

for i=1:length(S2_time)
    PO7_UR_epoch_values(:,i) = PO7(floor(UR_epoch_matrix(1,i)):floor(UR_epoch_matrix(3,i)))-UR_baseline(1,i);
end

for i=1:length(S1_time)
    PO8_UL_epoch_values(:,i) = PO8(floor(UL_epoch_matrix(1,i)):floor(UL_epoch_matrix(3,i)))-UL_baseline(2,i);
end

for i=1:length(S2_time)
    PO8_UR_epoch_values(:,i) = PO8(floor(UR_epoch_matrix(1,i)):floor(UR_epoch_matrix(3,i)))-UR_baseline(2,i);
end
    

%% Calculate the template (or mean shape) for both tasks in both channels : you should obtain 4 templates in total!

PO7_UL_Template = sum(PO7_UL_epoch_values,2)/length(S1_time);
PO7_UR_Template = sum(PO7_UR_epoch_values,2)/length(S2_time);
PO8_UL_Template = sum(PO8_UL_epoch_values,2)/length(S1_time);
PO8_UR_Template = sum(PO8_UR_epoch_values,2)/length(S2_time);

%% Calculate the overal SNR for the UL and UR tasks in both channels. Show the results in a table

% noise power 

noise_PO7_UL = 1/(length(PO7_UL_epoch_values)*(length(S1_time)-1))*sum(sum(((PO7_UL_epoch_values-PO7_UL_Template).^2),1),2);
noise_PO7_UR = 1/(length(PO7_UR_epoch_values)*(length(S2_time)-1))*sum(sum(((PO7_UR_epoch_values-PO7_UR_Template).^2),1),2);
noise_PO8_UL = 1/(length(PO8_UL_epoch_values)*(length(S1_time)-1))*sum(sum(((PO8_UL_epoch_values-PO8_UL_Template).^2),1),2);
noise_PO8_UR = 1/(length(PO8_UR_epoch_values)*(length(S2_time)-1))*sum(sum(((PO8_UR_epoch_values-PO8_UR_Template).^2),1),2);


% signal power 

signal_PO7_UL = 1/(length(PO7_UL_epoch_values))*sum(((PO7_UL_Template).^2),1)-(noise_PO7_UL/length(PO7_UL_epoch_values));
signal_PO7_UR = 1/(length(PO7_UR_epoch_values))*sum(((PO7_UR_Template).^2),1)-(noise_PO7_UR/length(PO7_UR_epoch_values));
signal_PO8_UL = 1/(length(PO8_UL_epoch_values))*sum(((PO8_UL_Template).^2),1)-(noise_PO8_UL/length(PO8_UL_epoch_values));
signal_PO8_UR = 1/(length(PO8_UR_epoch_values))*sum(((PO8_UR_Template).^2),1)-(noise_PO8_UR/length(PO8_UR_epoch_values));

% SNR

SNR_PO7_UL= signal_PO7_UL/noise_PO7_UL;
SNR_PO7_UR= signal_PO7_UR/noise_PO7_UR;
SNR_PO8_UL= signal_PO8_UL/noise_PO8_UL;
SNR_PO8_UR= signal_PO8_UR/noise_PO8_UR;

% Eucledian distance
D_PO7_UL = 1/(length(S1_time))* sum(sqrt(sum(((PO7_UL_epoch_values-PO7_UL_Template).^2),1)));
D_PO7_UR = 1/(length(S2_time))* sum(sqrt(sum(((PO7_UR_epoch_values-PO7_UL_Template).^2),1)));
D_PO8_UL = 1/(length(S1_time))* sum(sqrt(sum(((PO8_UL_epoch_values-PO8_UL_Template).^2),1)));
D_PO8_UR = 1/(length(S2_time))* sum(sqrt(sum(((PO8_UR_epoch_values-PO8_UL_Template).^2),1)));



%% Plot the required figures (see assignment)

% Plot of templates
time_end= length(PO7_UR_epoch_values) *1/fs;
labels = {'PO7 S1(UL)','PO7 S2(UR)','PO8 S1(UL)','PO8 S2(UR)'}'; 
x= 0:1/fs:(time_end-1/fs);

figure('Name','Plot of templates','NumberTitle','off');

subplot(3,2,1);plot(x,PO7_UL_Template,'b')
title('Template of PO7 S1(UL)')
xlabel('Time (s)')
ylabel('Amplitude in \muV')
subplot(3,2,2);plot(x,PO7_UR_Template,'g')
title('Template of PO7 S2(UR)')
xlabel('Time (s)')
ylabel('Amplitude in \muV')
subplot(3,2,3);plot(x,PO8_UL_Template,'r')
title('Template of PO8 S1(UL)')
xlabel('Time (s)')
ylabel('Amplitude in \muV')
subplot(3,2,4);plot(x,PO8_UR_Template,'k')
title('Template of PO8 S2(UR)')
xlabel('Time (s)')
ylabel('Amplitude in \muV')


subplot(3,2,5);
hold on 
plot(x,PO7_UL_Template,'b')
plot(x,PO7_UR_Template,'g')
plot(x,PO8_UL_Template,'r')
plot(x,PO8_UR_Template,'k')
xlabel('Time (s)')
ylabel('Amplitude in \muV')
title('Plot of all ERPs')
legend(labels)
hold off 


% Plot of ERP's


figure('Name','Plot of ERPs','NumberTitle','off');

subplot(2,2,1);
hold on 
plot(x,PO7_UL_epoch_values,'c','LineWidth',0.01)
plot(x,PO7_UL_Template,'-b','LineWidth',3)
title('ERPs with template overlap for PO7 S1(UL)')
xlabel('Time (s)')
ylabel('Amplitude in \muV')
hold off

subplot(2,2,2);
hold on 
plot(x,PO7_UR_epoch_values,'c','LineWidth',0.01)
plot(x,PO7_UR_Template,'-g','LineWidth',3)
title('ERPs with template overlap for PO7 S2(UR)')
xlabel('Time (s)')
ylabel('Amplitude in \muV')
hold off

subplot(2,2,3);
hold on 
plot(x,PO8_UL_epoch_values,'c','LineWidth',0.01)
plot(x,PO8_UL_Template,'-r','LineWidth',3)
title('ERPs with template overlap for PO8 S1(UL)')
xlabel('Time (s)')
ylabel('Amplitude in \muV')
hold off

subplot(2,2,4);
hold on 
plot(x,PO8_UR_epoch_values,'c','LineWidth',0.01)
plot(x,PO8_UR_Template,'-k','LineWidth',3)
title('ERPs with template overlap for PO8 S2(UR)')
xlabel('Time (s)')
ylabel('Amplitude in \muV')
hold off

% Table plotting

SNR_values = [SNR_PO7_UL,SNR_PO7_UR,SNR_PO8_UL,SNR_PO8_UR];
D_values = [D_PO7_UL,D_PO7_UR,D_PO8_UL,D_PO8_UR];
data= zeros(2,4);
data(1,:)=SNR_values;
data(2,:)=D_values;
f=figure('Name','SNR And Euclidean table','NumberTitle','off');
uit = uitable(f,'Data',data','Position',[20 20 710 310],'ColumnName',{'SNR values','Euclidean distance'},'RowName',labels);