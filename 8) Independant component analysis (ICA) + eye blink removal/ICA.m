-------------------------

    clear; close all; clc; 
    addpath ./FastICA_25/ % make sure MATLAB finds the fastICA-algorithm
    rng('default');



    load('MixtureECG');

    % --- Define the time axis for plotting

    t = (0:size(S,2)-1)/fs; 

    % --- Normalize the sources: subtract the mean and divide by the standard
    % deviation

    S = zscore(S,0,2);

  
    figure('Name','Normalized source signals','NumberTitle','off');
    for m = 1:size(S,1)
        subplot(size(S,1),1,m);
        plot(t,S(m,:));
        title(['Source #',num2str(m)]);
        xlabel('Time [s]');
    end

    %% 1.1.1 Measuring four mixtures

    % --- Mix the source signals

    X = A*S;
   
    % --- Plot the mixtures (use the same style of figure as above)

    figure('Name','Mixtures','NumberTitle','off');
    for m = 1:size(X,1)
        subplot(size(X,1),1,m);
        plot(t,X(m,:));
        title(['Sensor #',num2str(m)]);
        xlabel('Time [s]');
    end

    % --- Demix the mixtures using fastICA

    [Y,B,W,mm]= fastica(X);

    % --- Normalize the estimated sources: subtract the mean and divide by the 
    % standard deviation
 
    Y = zscore(Y,0,2);
   
    % --- Plot the estimated sources (use the same style of figure as above)

    figure('Name','Esdtimated source (4 mixtures)','NumberTitle','off');
    for m = 1:size(Y,1)
        subplot(size(Y,1),1,m);
        plot(t,Y(m,:));
        title(['Estimated source #',num2str(m)]);
        xlabel('Time [s]');
    end
        
% --- Identification of the estimated sources: automatically match the
% estimated sources with the corresponding original sources by using
% correlation as a similarity metric.


match_mat = zeros(4,4);

% Store pairwise correlation between each source and estimate 

for j=1:length(match_mat) 
    for i=1:length(match_mat)
        match_mat(i,j) = (corr(S(j,:)',Y(i,:)'));
    end
end

% Find max index of each column (for plotting of corresponding
% source-estimate pair)

[~,iMaxCorr] = max(abs(match_mat));

% --- Plot the matched pairs of estimated and original source next to each other.

figure('Name','Original vs Estimated source(4 mix)','NumberTitle','off');
for m = 1:size(S,1)
    subplot(size(S,1),2,2*m-1); 
    % plot Original source
    plot(t,S(m,:));
    xlabel('Time [s]'); title(['Source #' num2str(m)]);
    subplot(size(S,1),2,2*m); 
    % plot corresponding estimated source
    plot(t,Y(iMaxCorr(m),:))
    xlabel('Time [s]'); title(['Estimate #' num2str(iMaxCorr(m))]);
end

% --- Quality assessment of the estimated sources: calculate the
% root-mean-square error (RMSE) between each estimated source and its
% matched original source. Make sure the signals are standardized and that 
% they have the same sign (e.g., by picking the minimal RMSE for both 
% possible signs).

rmse = zeros(4,1);
Yn=-Y;

for i=1:length(rmse)
        a = sqrt(mean((S(i,:)-Y(iMaxCorr(i),:)).^2));
        b = sqrt(mean((S(i,:)-Yn(iMaxCorr(i),:)).^2));
        rmse(i,1) = min([a,b]);
end



% --- Artifact removal using ICA: remove the white Gaussian noise and
% sawtooth signal from the mixtures by modifying the mixing matrix.

sawtooth = S(2,:); gaussNoise = S(4,:); % the sawtooth signal and Gaussian noise source are saved in corresponding variables


 %Find corresponding components to noise from ICA output through correlation 
 denoise_mat = zeros(4,2);

for i=1:length(denoise_mat) 
    
      denoise_mat(i,1) = corr(Y(i,:)',sawtooth');
      denoise_mat(i,2) = corr(Y(i,:)',gaussNoise');
   
end

% Find index corresponding to the row component of the output mixing matrix
% of fastICA that needs to be zeros count

[~,idx_noise] = max(abs(denoise_mat));

%% 
% zero out noise components 

B(:,[idx_noise(1),idx_noise(2)])=0;


% Use manipulated mixing matrix to denoise the signal 
 X_denoise= B*Y;

% TIP: x(t) = Ay(t). Should you normalize the estimated sources in y(t)?

 %X_denoise = zscore( X_denoise,0,2);
% --- Plot the original mixtures next to the cleaned mixtures (use the same style of figure as above).

figure('Name','Originapl vs Estimated denoised source (4 mix)','NumberTitle','off');
for m = 1:size( X,1)
    subplot(size( X,1),2,2*m-1); 
    % plot Original source
    plot(t, X(m,:));
    xlabel('Time [s]'); title(['Source #' num2str(m)]);
    subplot(size(S,1),2,2*m); 
    % plot corresponding estimated source
    plot(t, X_denoise(m,:))
    xlabel('Time [s]'); title(['Estimate #' num2str(iMaxCorr(m))]);
end

    
%% 1.1.2 Measuring 30 mixtures

% --- Mix the source signals, without additional noise

X2 = A2*S;

% --- fastICA without additional noise: apply fastICA to mixture X2,
% normalize the estimated sources and match them. Compute the eigenvalues
% of the covariance matrix to support your analysis.


% ICA
[Y2,B2,W2]= fastica(X2); 

% normalization
 Y2 = zscore(Y2,0,2); 

% matching and plotting
match_mat2 = zeros(4,4);

for j=1:length(match_mat2) 
    for i=1:length(match_mat2)
        match_mat2(i,j) = (corr(S(j,:)',Y2(i,:)'));
    end
end

[~,iMaxCorr2] = max(abs(match_mat2));

% --- Plot the matched pairs of estimated and original source next to each other.

figure('Name','Original vs Estimated source (30 mix)','NumberTitle','off');
for m = 1:size(S,1)
    subplot(size(S,1),2,2*m-1); 
    % plot Original source
    plot(t,S(m,:));
    xlabel('Time [s]'); title(['Source #' num2str(m)]);
    subplot(size(S,1),2,2*m); 
    % plot corresponding estimated source
    plot(t,Y2(iMaxCorr2(m),:))
    xlabel('Time [s]'); title(['Estimate #' num2str(iMaxCorr2(m))]);
end

% RMSE computation

rmse2 = zeros(4,1);
Yn2=-Y2;

for i=1:length(rmse2)

        a = sqrt(mean((S(i,:)-Y2(iMaxCorr2(i),:)).^2));
        b = sqrt(mean((S(i,:)-Yn2(iMaxCorr2(i),:)).^2));
        rmse2(i,1) = min([a,b]);
        
end 

% eigenvalue computation

a2 = cov(X2');
[V2,D2] = eig(a2);
eig_diag2 = diag(D2);
eig_val2 = eig_diag2(diag(D2)>1);
%Result: 4 non zero eigenvalues  corresponding to the number of components detected by fastICA 
[E2, D2_ica] = fastica(X2, 'only', 'pca');
round(diag(D2_ica))== round(eig_val2); %Confirmation that eigenvalues are indeed equal 

% --- Mix the source signals, now with additional noise

X3 = A2*S + 0.75*(rand(size(A2,1),length(S))-0.5); % note: -0.5 to center the noise

% --- fastICA with additional noise: apply fastICA to mixture X3,
% normalize the estimated sources and match them. Compute the eigenvalues
% of the covariance matrix to support your analysis.

% TIP: check the preprocessing and the output of the fastICA algorithm to 


% ICA
a3 = cov(X3');
[V3,D3] = eig(a3);
eig_diag3 = diag(D3);

[Y3,B3,W3]= fastica(X3);

%%

% matching and plotting
match_mat3 = zeros(30,4);
[ro,co] = size(match_mat3);

for j=1:co
    for i=1:ro
        match_mat3(i,j) = (corr(S(j,:)',Y3(i,:)'));
    end
end

[~,iMaxCorr3] = max(abs(match_mat3));

% --- Plot the matched pairs of estimated and original source next to each other.

figure('Name','Original vs Estimated source (30 mix+ random noise)','NumberTitle','off');
for m = 1:size(S,1)
    subplot(size(S,1),2,2*m-1); 
    % plot Original source
    plot(t,S(m,:));
    xlabel('Time [s]'); title(['Source #' num2str(m)]);
    subplot(size(S,1),2,2*m); 
    % plot corresponding estimated source
    plot(t,Y3(iMaxCorr3(m),:))
    xlabel('Time [s]'); title(['Estimate #' num2str(iMaxCorr3(m))]);
end


% RMSE computation
rmse3 = zeros(4,1);
Yn3=-Y3;

for i=1:length(rmse3)

        a = sqrt(mean((Y3(iMaxCorr3(i),:)-S(i,:)).^2));
        b = sqrt(mean((S(i,:)-Yn3(iMaxCorr3(i),:)).^2));
        rmse3(i,1) = min([a,b]);
        
end 

% --- Dimensionality reduction as a preprocessing step for ICA: improve on
% the previous source separation results by reducing the dimensionality
% using PCA.

% compute and plot the eigenvalues of the covariance matrix 
figure('Name','eigenvalues of the covariance matrix ','NumberTitle','off');
plot([1:length(eig_diag3)],flip(eig_diag3),'-o','MarkerFaceColor','r','MarkerEdgeColor','k')
ylabel('Eigenvalues')
xlabel('Subscript of Eigenvalues')


    
% reduce the dimensionality of X3 accordingly, using PCA.

[coeff,score,latent,~,~, mu]= pca(X3');  %Latent gives eigendecomposition of covariance matrix and coeff contain eigenvectros but in reverse order 

coeff = coeff(:,[1:4]);  % Use first 4 eigenvalues as inferred from the plotted graph  
X3_pca =coeff'*X3;   %Project data (Reduce dimensions) onto the set of the first 4 principal components 
X_reconstruct = X3_pca'*coeff'; %reconstruct the signal back to the sensor space 
X_reconstruct=X_reconstruct';

% reapply ICA
[Y3_pca,~,~]= fastica(X_reconstruct);

% normalization
Y3_pca = zscore(Y3_pca,0,2); 

% matching and plotting
match_mat3_pca = zeros(4,4);

for j=1:length(match_mat3_pca) 
    for i=1:length(match_mat3_pca)
        match_mat3_pca(i,j) = (corr(S(j,:)',Y3_pca(i,:)'));
    end
end

[~,iMaxCorr3_pca] = max(abs(match_mat3_pca));

figure('Name','Original vs Estimated source_reduced dimensions (30 mix+ random noise)','NumberTitle','off');
for m = 1:size(S,1)
    subplot(size(S,1),2,2*m-1); 
    % plot Original source
    plot(t,S(m,:));
    xlabel('Time [s]'); title(['Source #' num2str(m)]);
    subplot(size(S,1),2,2*m); 
    % plot corresponding estimated source
    plot(t,Y3_pca(iMaxCorr3_pca(m),:))
    xlabel('Time [s]'); title(['Estimate #' num2str(iMaxCorr3_pca(m))]);
end

% RMSE computation
rmse3_pca = zeros(4,1);
Yn3_pca=-Y3_pca;

for i=1:length(rmse3_pca)

        a = sqrt(mean((S(i,:)-Y3_pca(iMaxCorr3_pca(i),:)).^2));
        b = sqrt(mean((S(i,:)-Yn3_pca(iMaxCorr3_pca(i),:)).^2));
        rmse3_pca(i,1) = min([a,b]);
        
end 
