function [features, labels] = main_sleep(filename)
% Perform feature extraction & create sleep labels 

% INPUT:
% filename: name of file containing data

% OUTPUT: 
% features: feature matrix (number of segments x number of features)
% labels: sleep label for each 30s EEG segment

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Load the data
load(filename);

% --- Normalize the EEG signal: remove the DC component and divide by the
% standard deviation
data = zscore(EEG.data); 

% --- Segment the signal into nonoverlapping windows of 30s
len_seg = 30*EEG.fs;
nsegments = length(data)/(len_seg);
mat = zeros(nsegments,len_seg);
label_mat=zeros(nsegments,1);

j=1;
for i=1:nsegments
    mat(i,:)=data(j:len_seg*i);
    label_mat(i,:) = mode(EEG.labels(j:len_seg*i));
    j=len_seg*i+1;
end 

% --- Estimate PSD and extract relative power in delta, theta, alpha & beta
% band

[pxx,f] = pwelch(mat',1500,0,[],EEG.fs,'onesided');

x=[0.5 4 8 13 35]; %Frequency band
freq_ind = round(length(f)*x/(EEG.fs/2))+1; %Frequency indices 
mat2 = zeros(nsegments,4);

for i =1:nsegments
    c= pxx(:,i);
    total_power = sum(c(freq_ind(1):freq_ind(end)));
    delta = sum(c(freq_ind(1):freq_ind(2)))/total_power;
    theta = sum(c(freq_ind(2)+1:freq_ind(3)))/total_power;
    alpha = sum(c(freq_ind(3)+1:freq_ind(4)))/total_power;
    beta = sum(c(freq_ind(4)+1:freq_ind(5)))/total_power;

    mat2(i,1) = delta;
    mat2(i,2) = theta;
    mat2(i,3) = alpha;
    mat2(i,4) = beta;
end

% --- Construct sleep label for each 30s EEG segment

labels = label_mat; 

% --- Perform feature normalization (if needed)

mat2 = zscore(mat2,0,1); 

features = mat2; 
end


