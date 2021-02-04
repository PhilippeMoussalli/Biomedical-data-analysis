function [label_test] = knn(train,label_train,test,k)
% Perform kNN classification

% Input: - train: features of the training set
%        - label_train: labels of the training set
%        - test: feature of the test sample
%        - k: number of nearest neighbours to take into account

% Output: label_test: label of the test sample

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Compute distance to all other samples

dist_mat = zeros(1,length(train));

for i=1:length(train)
    d1 = sqrt((test-train(i,:)) *(test-train(i,:))');
    dist_mat(1,i)= d1;
end

[~,idx_asc] = sort(dist_mat,'ascend');

% --- Assign the label of majority of k closest samples

label_test = mode(label_train(idx_asc(1:k)));

end

