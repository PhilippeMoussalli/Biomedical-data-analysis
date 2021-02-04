-----------------------------------

%% Initialization

clear; close all; clc; 
dwtmode('per','nodisp'); % periodic extension for the wavelet decomposition
load('Wavelet_1_data.mat');

%% The wavelet decomposition of the ECG signal

% --- Use the Daubechies wavelet of order four (db4) with five levels to 
% decompose the ECG signal

[c_ecg,l_ecg] = wavedec(ecg,5,'db4'); %Wavelet Decomposition 
%rec_ecg = waverec(c_ecg,l_ecg,'db4'); %Wavelet Reconstruction 

A5 = appcoef(c_ecg,l_ecg,'db4',5); %Approximation Coefficient Level 5

% --- Plot the original ECG signal, together with the approximation and
% detail coefficients, in one subplot (stacked below each other) in the
% following order (top to bottom): (ECG,A5,D5,D4,D3,D2,D1).

figure(1)
subplot(7,1,1);
% plot ECG signal
plot(ecg)
hold on;
% plot reconstructed signal
%plot(rec_ecg)
title('ECG signal'); axis tight;
xlabel('Samples');
ylabel('a.u.');

subplot(7,1,2);
% plot approximation coefficients using stem
stem(A5) 
title('A5'); axis tight; 
xlabel('Samples');
ylabel('a.u.');

j=5;
for i = 1:5
    subplot(7,1,i+2);
    % plot detail coefficients using stem and give it the appropriate title
    D = detcoef(c_ecg,l_ecg,j);
    stem(D); 
    title(['D',sprintf('%d',j)]); axis tight;
    xlabel('Samples');
    ylabel('a.u.'); 
    j = j-1; 
end


%% Compression of the ECG signal using 10, 25, 50 and 100 coefficients

% --- Keep only M = 10,25,50,100 wavelet coefficients with the largest
% absolute value to compress the signal. Plot the results as above and
% compute for each M the corresponding RMSE and CR.

% TIP: code in a modular way, e.g., by using loops and functions (use the
% same plotting code as above)

M_values = [10 25 50 100]; 
RMSE_mat = zeros(4,1);
CR_mat = zeros(4,1);

for i = 1:4
[recons,coeff_total] = reconsWavelet(c_ecg,l_ecg,M_values(i));
RMSE_mat(i,1) = sqrt(mean((c_ecg - coeff_total).^2));
CR_mat(i,1) = length(ecg) / M_values(i);

figure(i+1)
subplot(7,1,1);
% plot ECG signal
plot(ecg,'b')
hold on;
% plot reconstructed signal
plot(recons,'r')
title('ECG signal'); axis tight;
xlabel('Samples');
ylabel('a.u.');
legend('Original Signal','Reconstructed Signal');

subplot(7,1,2);
% plot approximation coefficients using stem
stem(coeff_total(1:l_ecg(1))) 
title('A5'); axis tight;
xlabel('Samples');
ylabel('a.u.');

j=1;
k=5;
    for p = 1:5
        subplot(7,1,p+2);
        % plot detail coefficients using stem and give it the appropriate title
        stem(coeff_total(l_ecg(j+1)+1:l_ecg(j+2)+1)); 
        title(['D',sprintf('%d',k)]); axis tight;
        xlabel('Samples');
        ylabel('a.u.'); 
        j = j+1;
        k = k-1;
    end

end 

compression_values = {'10';'25';'50';'100'};
column_name = {'Compression Sizes','RMSE','CR'};
T_compression = table(compression_values,RMSE_mat(:,1),CR_mat(:,1));
T_compression.Properties.VariableNames = column_name;
T_compression

%%Additional Figure for Report 

[recons,coeff_total] = reconsWavelet(c_ecg,l_ecg,10);

figure
% plot ECG signal
plot(ecg,'b')
hold on;
% plot reconstructed signal
plot(recons,'r')
title('ECG signal / 10 coefficients compression'); axis tight;
xlabel('Samples');
ylabel('a.u.');
legend('Original Signal','Reconstructed Signal');

