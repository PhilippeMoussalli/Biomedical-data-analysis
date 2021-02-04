%--------------------------------------------------------------------------
%
% EX 4.2 : Unsupervised classification
%
%--------------------------------------------------------------------------

clear all;
close all;
clc;

% --- Extract features and labels of the EEG data
filename = 'EEG_sleep.mat';
[data, labels] = main_sleep(filename);
features = {'Relative Delta power','Relative Theta power','Relative Alpha power','Relative Beta power'};
sleep_stages = {'N3','N2','N1','REM','Wake'};
[feature_mat, labels] = main_sleep(filename);
%% --- Perform k-means clustering with 5 clusters and visualize the result

% Perform 5-means clustering
num_clusters = 5;
cmb = combnk(1:4,2);
cent_mat = zeros(num_clusters*length(cmb),2);
idx_mat = zeros(length(data),length(cmb));
sumd_mat = zeros(num_clusters,length(cmb));
j=1;

for i=1:length(cmb)
    X= [feature_mat(:,cmb(i,1)),feature_mat(:,cmb(i,2))];
    [idx,C,sumd,D] = kmeans(X,num_clusters);
    idx_mat(:,i) = idx;
    cent_mat((j:j+num_clusters-1),:)= C;
    j= (j+num_clusters);
    sumd_mat(:,i) = sumd;
end

% Plot result of clustering
j=1;
for i=1:length(cmb)
    subplot(3,2,i)
    hold on 
    gscatter(data(:,cmb(i,1)),data(:,cmb(i,2)),idx_mat(:,i),'rbmgy')
    xlabel(features(cmb(i,1)))
    ylabel(features(cmb(i,2)))
    plot(cent_mat(j:j+num_clusters-1,1),cent_mat(j:j+num_clusters-1,2),'kx','LineWidth',2);
    j= (j+num_clusters);
    legend('1','2','3','4','5','centroid')
    xlim([-4 4])
end
%% --- Perform k-means clustering with 3 clusters and visualize the result

% Perform 3-means clustering
num_clusters_2 = 3;
cent_mat_2 = zeros(num_clusters_2*length(cmb),2);
idx_mat_2 = zeros(length(data),length(cmb));
sumd_mat_2 = zeros(num_clusters_2,length(cmb));
j=1;

for i=1:length(cmb)
    X= [feature_mat(:,cmb(i,1)),feature_mat(:,cmb(i,2))];
    [idx_2,C_2,sumd_2] = kmeans(X,num_clusters_2);
    idx_mat_2(:,i) = idx_2;
    cent_mat_2((j:j+num_clusters_2-1),:)= C_2;
    j= (j+num_clusters_2);
    sumd_mat_2(:,i) = sumd_2;
end

% Plot result of clustering

figure
j=1;
for i=1:length(cmb)
    subplot(3,2,i)
    hold on 
    gscatter(data(:,cmb(i,1)),data(:,cmb(i,2)),idx_mat_2(:,i),'rbm')
    xlabel(features(cmb(i,1)))
    ylabel(features(cmb(i,2)))
    plot(cent_mat_2(j:j+num_clusters_2-1,1),cent_mat_2(j:j+num_clusters_2-1,2),'kx','LineWidth',2);
    j= (j+num_clusters_2);
    legend('1','2','3','Centroid')
    xlim([-4 4])
end


% Table plotting for distance between center and clusters (K=5)
column_name = {'Cluster number','Beta/Alpha','Beta/Theta','Alpha/Theta','Beta/Delta','Alpha/Delta','Theta/Delta'};
distance_t_5 = table({1,2,3,4,5}',sumd_mat(:,1),sumd_mat(:,2),sumd_mat(:,3),sumd_mat(:,4),sumd_mat(:,5),sumd_mat(:,6));
distance_t_5.Properties.VariableNames = column_name;
distance_t_5

% Try out different methods for intialization of centroid 


cent_mat_3 = zeros(num_clusters*length(cmb),2);
idx_mat_3 = zeros(length(data),length(cmb));
sumd_mat_3 = zeros(num_clusters,length(cmb));
j=1;

for i=1:length(cmb)
    X= [feature_mat(:,cmb(i,1)),feature_mat(:,cmb(i,2))];
    [idx_3,C_3,sumd_3] = kmeans(X,num_clusters,'Start','sample');
    idx_mat_3(:,i) = idx_3;
    cent_mat_3((j:j+num_clusters-1),:)= C_3;
    j= (j+num_clusters);
    sumd_mat_3(:,i) = sumd_3;
end

% Plot result of clustering
figure
j=1;
for i=1:length(cmb)
    subplot(3,2,i)
    hold on 
    gscatter(data(:,cmb(i,1)),data(:,cmb(i,2)),idx_mat_3(:,i),'rbmgy')
    xlabel(features(cmb(i,1)))
    ylabel(features(cmb(i,2)))
    plot(cent_mat_3(j:j+num_clusters-1,1),cent_mat_3(j:j+num_clusters-1,2),'kx','LineWidth',2);
    j= (j+num_clusters);
    legend('1','2','3','4','5','centroid')
    xlim([-4 4])
end
