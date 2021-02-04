
clear all;
close all;
clc;

% --- Extract features & labels of EEG data
features = {'Relative Delta power','Relative Theta power','Relative Alpha power','Relative Beta power'};
sleep_stages = {'N3','N2','N1','REM','Wake'};

% --- Load the features
load('EEG_sleep_supervised.mat');

% --- Perform feature normalization (if needed)
 feat_test  = zscore(feat_test);
 feat_train = zscore(feat_train);


%% Part 1: on training data

% --- Perform leave-one-out cross-validation and predict the label of the
% test sample using the distance function, kNN and LDA
pred_label_df = zeros(1,length(labels_train));
pred_label_k1 = zeros(1,length(labels_train));
pred_label_k3 = zeros(1,length(labels_train));
pred_label_k5 = zeros(1,length(labels_train));
pred_label_k7 = zeros(1,length(labels_train));
pred_label_LDA = zeros(1,length(labels_train));

%Data partitoning

for i=1:length(feat_train)
    test_LOO = feat_train(i,:);
    train_LOO = feat_train;
    labels_train_LOO = labels_train;
    train_LOO(i,:)=[];
    labels_train_LOO(i)=[];

   %%% --- Distance function

   % complete and apply the function 'df.m'
   
    [label_test_df] = df(train_LOO,labels_train_LOO,test_LOO);
    pred_label_df(1,i) = label_test_df;
   
   %%% --- KNN

   % complete and apply the function 'knn.m'
   [label_test_k1] = knn(train_LOO,labels_train_LOO,test_LOO,1);
   [label_test_k3] = knn(train_LOO,labels_train_LOO,test_LOO,3);
   [label_test_k5] = knn(train_LOO,labels_train_LOO,test_LOO,5);
   [label_test_k7] = knn(train_LOO,labels_train_LOO,test_LOO,7);
   
   pred_label_k1(1,i) = label_test_k1;
   pred_label_k3(1,i) = label_test_k3;
   pred_label_k5(1,i) = label_test_k5;
   pred_label_k7(1,i) = label_test_k7;
   
   %%% --- LDA

   Mdl = fitcdiscr(train_LOO,labels_train_LOO);
   label_test_LDA = predict(Mdl,test_LOO);
   pred_label_LDA(1,i) = label_test_LDA;
   
end


% --- Compute performance (accuracy, sensitivity and specificity)

[acc_df,sens_df,spec_df] = performance(pred_label_df,labels_train);
[acc_k1,sens_k1,spec_k1] = performance(pred_label_k1,labels_train);
[acc_k3,sens_k3,spec_k3] = performance(pred_label_k3,labels_train);
[acc_k5,sens_k5,spec_k5] = performance(pred_label_k5,labels_train);
[acc_k7,sens_k7,spec_k7] = performance(pred_label_k7,labels_train);
[acc_LDA,sens_LDA,spec_LDA] = performance(pred_label_LDA,labels_train);


% Make table showing the performance of the different methods

Classification_method = {'Distance function';'KNN_1';'KNN_3';'KNN_5';'KNN_7';'LDA'};
accuracy = [acc_df;acc_k1;acc_k3;acc_k5;acc_k7;acc_LDA];
sensitivity = [sens_df;sens_k1;sens_k3;sens_k5;sens_k7;sens_LDA];
specificity = [spec_df;spec_k1;spec_k3;spec_k5;spec_k7;spec_LDA];
classifier_T= table(Classification_method,accuracy,sensitivity,specificity);
classifier_T


%% Part 2: on new, unseen test data

% --- Train LDA classifier using training data and apply on test data

% --- Make ROC curve & compute area under curve

Mdl2 = fitcdiscr(feat_train,labels_train); 
[label_test_LDA2,score] = predict(Mdl2,feat_test); 
roc_mat = [score(:,2) labels_test'];

[acc_df,sens_df,spec_df] = performance(label_test_LDA2,labels_test);

spacing = linspace(0,1,1000);
roc_plot_mat = zeros(length(spacing)-1,2);
roc_plot_mat2 = zeros(length(spacing)-1,2);


%Compute ROC
for i=1:length(spacing)
    idx_thresh= find(roc_mat(:,1)>spacing(i));
    test_label_roc = zeros(1,length(roc_mat));
    test_label_roc(idx_thresh)=1;
    [~,sens,spec] = performance(test_label_roc,roc_mat(:,2));
    roc_plot_mat(i,1)=sens;
    roc_plot_mat(i,2)=1-spec;
end


% Compute and display AUC
AUC= abs(trapz(roc_plot_mat(:,2),roc_plot_mat(:,1)));

% --- Mark requested working point on ROC curve
ind_wp=find(round(roc_plot_mat(:,2),2,'significant')==0.10);
sens_wp = roc_plot_mat(ind_wp,1);

figure
hold on 
area(roc_plot_mat(:,2),roc_plot_mat(:,1))
plot([0,1],[0,1],'--k')
plot(roc_plot_mat(ind_wp,2),roc_plot_mat(ind_wp,1),'rx','LineWidth',5)
legend('ROC','Random classifier','Working point')
xlabel('FPR')
ylabel('TPR')
title(strcat('ROC Curve of LDA classifier on test data. AUC =',sprintf('%.2f',AUC)))
