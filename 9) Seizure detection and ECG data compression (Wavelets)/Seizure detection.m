
%% Initialization

clear; close all; clc; 


load('Seizure_data.mat');

%% Data exploration

% --- Plot all non-seizure segments in one subplot, with all seizure
% segments in the subplot below.

t = (1/fs:1/fs:(length(dataTrain)/fs));
non_seizure = dataTrain(find(labelsTrain==0),:);
seizure = dataTrain(find(labelsTrain==1),:);

figure(1) %Low vs. High frequency observed 
subplot(2,1,1);
plot(t,non_seizure); 
xlabel('Time [s]'); ylabel('a.u.'); title('Non-seizure'); ylim([-800,800]);
subplot(2,1,2); 
plot(t,seizure); 
xlabel('Time [s]'); ylabel('a.u.'); title('Seizure'); ylim([-800,800]);


% --- Do you think we should standardize the raw data? If so, add your code
% here.

data_std = std(dataTrain,0,2); 
data_mean = mean(dataTrain,2);

dataTrain = dataTrain - data_mean;

%% Full signal energy as a feature

% --- Calculate the signal energy for all training segments.

full_energy = sum(dataTrain.^2,2);

%% Wavelet band energies as features

% --- Decompose each training segment with the wavelet decomposition, using
% the db4 wavelet and calculate the energy contained in the coefficients of
% each wavelet band.

wave_energy = zeros(size(dataTrain,1),6);

for i=1:size(dataTrain,1)
[c_eeg,l_eeg] = wavedec(dataTrain(i,:),5,'db4'); 
A5 = appcoef(c_eeg,l_eeg,'db4',5);                      %Approximation coefficients
wave_energy(i,1) = sum(A5.^2);                          %Approximation coefficients energy
    
    j=5;
    for k=1:5                                           %Detail coefficients and their energy
        D = detcoef(c_eeg,l_eeg,j); 
        wave_energy(i,k+1) = sum(D.^2);
        j=j-1;
    end 

end

%% Visualization of the features

% --- Visualize all features (full signal energy + all wavelet band 
% energies) using boxplots (7 subplots in total).

figure;
subplot(1,7,1); 
no_seiz_full = full_energy(find(labelsTrain==0),1);
seiz_full = full_energy(find(labelsTrain==1),1);
x_full=[no_seiz_full;seiz_full];
g0_full = zeros(length(no_seiz_full),1);
g1_full = ones(length(seiz_full),1);
g_full = [g0_full;g1_full];
boxplot(x_full,g_full)
set(gca, 'XTickLabel', {'Non-Seizure','Seizure'});
ylabel('Energy (a.u.)');
title('Full energy'); 
for j = 1:6
    subplot(1,7,j+1); 
    no_seiz_wave = wave_energy(find(labelsTrain==0),j);
    seiz_wave = wave_energy(find(labelsTrain==1),j);
    x_wave = [no_seiz_wave;seiz_wave];
    g0_wave = zeros(length(no_seiz_wave),1);
    g1_wave = ones(length(seiz_wave),1);
    g_wave = [g0_wave;g1_wave];
    boxplot(x_wave,g_wave)
    set(gca, 'XTickLabel', {'Non-Seizure','Seizure'});
    % give each subplot a title
    titles = {'A5','D5','D4','D3','D2','D1'};
    title(sprintf('%s',titles{j}));
end


%% Classification of epileptic seizure activity

% --- Classification with full energy as only features

% train a classifier using fitcdiscr using the full energy as only feature

Mdl_full = fitcdiscr(full_energy,labelsTrain);

% compute the full energy for the test set (similarly as before)

data_mean2 = mean(dataTest,2);    %Normalization 
dataTest = dataTest - data_mean2; 

full_energy2 = sum(dataTest.^2,2);

% predict the seizure activity using predict

[label_test_full,score_full] = predict(Mdl_full,full_energy2);

% --- Classification with wavelet subband energies as features

% train a classifier using fitcdiscr using the subband energies

Mdl_wave = fitcdiscr(wave_energy,labelsTrain);

% compute the subband energies for the test set (similarly as before)

wave_energy2 = zeros(size(dataTest,1),6);

for i=1:size(dataTest,1)
[c_eeg,l_eeg] = wavedec(dataTest(i,:),5,'db4'); 
A5 = appcoef(c_eeg,l_eeg,'db4',5);                      %Approximation coefficients
wave_energy2(i,1) = sum(A5.^2);                          %Approximation coefficients energy
    
    j=5;
    for k=1:5                                           %Detail coefficients and their energy
        D = detcoef(c_eeg,l_eeg,j); 
        wave_energy2(i,k+1) = sum(D.^2);
        j=j-1;
    end 

end

% predict the seizure activity using predict

[label_test_wave,score_wave] = predict(Mdl_wave,wave_energy2);

%% Performance evaluation

% --- Evaluate your classifiers (accuracy, sensitivity, specificity)

[acc_full,sens_full,spec_full] = performance(label_test_full,labelsTest);
[acc_wave,sens_wave,spec_wave] = performance(label_test_wave,labelsTest);

    %Table to show results
Classification_method = {'Full Energy';'Subband Energies'};    
Accuracy = [acc_full;acc_wave];
Sensitivity = [sens_full;sens_wave];
Specificity = [spec_full;spec_wave];
classifier_T= table(Classification_method,Accuracy,Sensitivity,Specificity);
classifier_T

    %ROC curve to compare results

roc_mat_full = [score_full(:,2) labelsTest];
roc_mat_wave = [score_wave(:,2) labelsTest];

spacing = linspace(0,1,1000);
roc_plot_mat_full = zeros(length(spacing)-1,2);
roc_plot_mat_wave = zeros(length(spacing)-1,2);

%Compute ROC
for i=1:length(spacing)
    idx_thresh1= find(roc_mat_full(:,1)>spacing(i));                    %For full energy
    test_label_roc1 = zeros(1,length(roc_mat_full));
    test_label_roc1(idx_thresh1)=1;
    [~,sens1,spec1] = performance(test_label_roc1,roc_mat_full(:,2));
    roc_plot_mat_full(i,1)=sens1;
    roc_plot_mat_full(i,2)=1-spec1;

    idx_thresh2= find(roc_mat_wave(:,1)>spacing(i));                    %For subband energies 
    test_label_roc2 = zeros(1,length(roc_mat_wave));
    test_label_roc2(idx_thresh2)=1;
    [~,sens2,spec2] = performance(test_label_roc2,roc_mat_wave(:,2));
    roc_plot_mat_wave(i,1)=sens2;
    roc_plot_mat_wave(i,2)=1-spec2;
end

AUC1= abs(trapz(roc_plot_mat_full(:,2),roc_plot_mat_full(:,1)));
AUC2= abs(trapz(roc_plot_mat_wave(:,2),roc_plot_mat_wave(:,1)));

figure
hold on 
plot([0,1],[0,1],'--k')
plot(roc_plot_mat_full(:,2),roc_plot_mat_full(:,1),'b')
plot(roc_plot_mat_wave(:,2),roc_plot_mat_wave(:,1),'r')
legend('Random Classifier',sprintf('Full Energy. AUC = %.2f',AUC1),sprintf('Subband Energies. AUC = %.2f',AUC2))
xlabel('FPR')
ylabel('TPR')
title('ROC Curves of Different Classifiers');