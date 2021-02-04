function [label_test] = df(train,labels_train,test)
% Perform classification based on the distance function

% Input: - train: features of the training set
%        - label_train: labels of the training set
%        - test: feature of the test sample

% Output: label_test: label of the test sample

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Compute the prototypes

cent_mat = zeros(2,4);

j=0;
for i=1:2
    cent_mat(i,1) = mean(train(find(labels_train==j),1));
    cent_mat(i,2) = mean(train(find(labels_train==j),2));
    cent_mat(i,3) = mean(train(find(labels_train==j),3));
    cent_mat(i,4) = mean(train(find(labels_train==j),4));
    j=j+1;
end

% --- Compute the distance to both prototypes
% --- Assign label of the closest prototype

%Calulculate Euclidean distance from each class protoype 
d1 = sqrt((test-cent_mat(1,:)) *(test-cent_mat(1,:))');
d2 = sqrt((test-cent_mat(2,:)) *(test-cent_mat(2,:))');

     if d1>d2
       label_test= 1;
     else
      label_test= 0;
     end 
 
end



