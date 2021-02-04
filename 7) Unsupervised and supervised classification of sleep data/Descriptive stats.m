%--------------------------------------------------------------------------
%
% EX 4.1 : Descriptive statistics
%
%--------------------------------------------------------------------------

close all;
clear all;
clc;

% --- Extract features and labels of the EEG data
filename = 'EEG_sleep.mat';
load('EEG_sleep.mat'); 
time = (0:length(EEG.data)-1)/EEG.fs;
[feature_mat, labels] = main_sleep(filename);
features = {'Relative Delta power','Relative Theta power','Relative Alpha power','Relative Beta power'};
sleep_stages = {'N3','N2','N1','REM','Wake'};

%% --- Make a scatter plot of all possible feature combinations

cmb = combnk(1:4,2); % all possible feature combinations

% go with a for loop through all possible feature combinations and
% visualize the features (make 1 figure with 6 subplots)
figure
for i=1:length(cmb)
    subplot(3,2,i)
    hold on 
    scatter(feature_mat(find(labels==1),cmb(i,1)),feature_mat(find(labels==1),cmb(i,2)),'r');
    scatter(feature_mat(find(labels==2),cmb(i,1)),feature_mat(find(labels==2),cmb(i,2)),'b');
    scatter(feature_mat(find(labels==3),cmb(i,1)),feature_mat(find(labels==3),cmb(i,2)),'k');
    scatter(feature_mat(find(labels==4),cmb(i,1)),feature_mat(find(labels==4),cmb(i,2)),'g');
    scatter(feature_mat(find(labels==5),cmb(i,1)),feature_mat(find(labels==5),cmb(i,2)),'y');
    xlabel(features(cmb(i,1)))
    ylabel(features(cmb(i,2)))
    legend(sleep_stages)
end


%% --- Make boxplots for each feature
% (make 1 figure with 4 subplots)
figure
for i=1:4
    subplot(2,2,i)
    n3 = feature_mat(find(labels==1),i);
    n2 = feature_mat(find(labels==2),i);
    n1 = feature_mat(find(labels==3),i);
    REM = feature_mat(find(labels==4),i);
    wake = feature_mat(find(labels==5),i);
    x=[n3;n2;n1;REM;wake];
    g1 = repmat({'N3'},length(n3),1);
    g2 = repmat({'N2'},length(n2),1);
    g3 = repmat({'N1'},length(n1),1);
    g4 = repmat({'REM'},length(REM),1);
    g5 = repmat({'wake'},length(wake),1);
    g = [g1; g2; g3 ; g4 ; g5];
    boxplot(x,g)
    %ylim([-4 4])
    title(features(i))
end
    
%% --- Compute & display mean and standard deviation
% (Use a table with the features in the rows & mean and standard deviation
% for each class in the columns)
mat_mean=zeros(5,4);
mat_std=zeros(5,4);

for i=1:4
    n3 = feature_mat(find(labels==1),i);
    n2 = feature_mat(find(labels==2),i);
    n1 = feature_mat(find(labels==3),i);
    REM = feature_mat(find(labels==4),i);
    wake = feature_mat(find(labels==5),i);
    mean_arr=[mean(n3);mean(n2);mean(n1);mean(REM);mean(wake)];
    std_arr=[std(n3);std(n2);std(n1);std(REM);std(wake)];
    mat_mean(:,i) = mean_arr;
    mat_std(:,i) = std_arr;
end

column_name = {'Class','Relative Delta power','Relative Theta power','Relative Alpha power','Relative Beta power'};
T_mean = table(sleep_stages',mat_mean(:,1),mat_mean(:,2),mat_mean(:,3),mat_mean(:,4));
T_mean.Properties.VariableNames = column_name;
T_std = table(sleep_stages',mat_std(:,1),mat_std(:,2),mat_std(:,3),mat_std(:,4));
T_std.Properties.VariableNames = column_name;
T_mean
T_std