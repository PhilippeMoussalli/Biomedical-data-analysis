%--------------------------------------------------------------------------
%
% EX 2.1 : Frequency-domain and Wiener filtering of the ECG signal
%
%--------------------------------------------------------------------------

%% 1.0 Initialization
clearvars;
close('all');

% --- load the signal in _data_ex2_1.mat_
load('data_ex2_1.mat')
fs = ECG.fs;
data = ECG.data;

% --- normalize the data: subtract the mean and divide by the standard
% deviation
data = (data - mean(data))/std(data);
% ... 
% --- plot the normalized ECG signal (use a time axis with seconds as unit)
n = length(data);
t = [1/fs:1/fs:(n/fs)];
hz = linspace(0,fs/2,n/2);
%fig_handle = figure; % create a handle ('label') for the figure
figure
plot(data)

%% 1.1 Butterworth filtering

%--- specify filter settings

ord = [3,7];
fc = [30,40];


% --- filter the signal using the given settings for filter order and cut-off
% --- frequency (use the command _butter_)

filtered_butter = zeros (4,n);  % Matrix initialization

% write modular code!
counter = 1;
	for i=1:2  % order
        for j=1:2   % Cutoff frequency
        [b,a] = butter(ord(i),fc(j)/(fs/2));
        filtered_butter(counter,:)= filter(b,a,data);  
        counter = counter+1;
        [h,w]=freqz(b,a,floor(n/2));
        figure('Name','Frequency response of notch filter with zeroes only','NumberTitle','off');
        plot(hz(1:length(hz)),20*log10(abs(h)));
        xlabel('Frequency (Hz)')
        ylabel('Magnitude (dB)')
        title('frequency response plot')
        end 
    end
  
    

% --- plot the filtered output signals (overlaid on the original signal). 
% --- (use the code below to arrange the graphs)

relspacing = 0.3;
nrows = numel(ord); vspace = relspacing/(nrows+1); height = (1-(nrows+1)*vspace)/nrows;
relspacing = 0.2;
ncols = numel(fc); hspace = relspacing/(ncols+1); width = (1-(ncols+1)*hspace)/ncols;
figure('name','Butterworth filtering');
ax = [];
idx = 1; ax(idx) = subplot('position',[ hspace+mod(idx-1,ncols)*(width+hspace) , vspace+(nrows-ceil(idx/ncols))*(height+vspace) , width , height ])
    % <plot the filtered signal using the lowest order and lowest cut-off frequency and label the axes>
    plot(t, filtered_butter(1,:))
    axis('tight')
    xlabel('Time (s)')
    ylabel('Voltage (mV)')
    title(sprintf('N = %d , f_c = %d Hz',ord(ceil(idx/ncols)),fc(mod(idx-1,ncols)+1)),'fontsize',12)
    
idx = 2; ax(idx) = subplot('position',[ hspace+mod(idx-1,ncols)*(width+hspace) , vspace+(nrows-ceil(idx/ncols))*(height+vspace) , width , height ])
	% <plot the filtered signal using the lowest order and highest cut-off frequency and label the axes>
     plot(t, filtered_butter(2,:))
    axis('tight')
    xlabel('Time (s)')
    ylabel('Voltage (mV)')
    title(sprintf('N = %d , f_c = %d Hz',ord(ceil(idx/ncols)),fc(mod(idx-1,ncols)+1)),'fontsize',12)
    
idx = 3; ax(idx) = subplot('position',[ hspace+mod(idx-1,ncols)*(width+hspace) , vspace+(nrows-ceil(idx/ncols))*(height+vspace) , width , height ])
    % <plot the filtered signal using the highest order and lowest cut-off frequency and label the axes>
    plot(t, filtered_butter(3,:))
    axis('tight')
    xlabel('Time (s)')
    ylabel('Voltage (mV)')
    title(sprintf('N = %d , f_c = %d Hz',ord(ceil(idx/ncols)),fc(mod(idx-1,ncols)+1)),'fontsize',12)
    
idx = 4; ax(idx) = subplot('position',[ hspace+mod(idx-1,ncols)*(width+hspace) , vspace+(nrows-ceil(idx/ncols))*(height+vspace) , width , height ])
	% <plot the filtered signal using the highest order and highest cut-off frequency and label the axes>
    plot(t, filtered_butter(4,:))
    axis('tight')
    xlabel('Time (s)')
    ylabel('Voltage (mV)')
    title(sprintf('N = %d , f_c = %d Hz',ord(ceil(idx/ncols)),fc(mod(idx-1,ncols)+1)),'fontsize',12)
    
annotation('textarrow',[0.5-0.1 0.5+0.1],0.2*vspace+repmat(0.5*vspace+(nrows)*(height+vspace),1,2),'String','f_c','fontsize',20)
annotation('textarrow',repmat(0.2*hspace,1,2),[0.5+0.1 0.5-0.1],'String','N','fontsize',20)

linkaxes(ax); 

%% 1.2 Wiener filtering

% --- Construct a model of the desired signal

% Select a few points on the ECG signal to help you create a piecewise 
% linear model
fprintf('Go to the figure of the ECG signal and\nzoom in on a relevant portion of the signal\nin which you can see the desired waveforms.\n')
%figure(fig_handle) % this will bring you back to the figure with the ECG signal
%keyboard % this will pause the execution, giving you time to execute the tip below:
% TIP: First, zoom in on a part of the ECG signal, in which you 
% can view a proper heartbeat. When you are
% ready, press F5 or click 'Run' to continue

% once you have created a 'good' model, save it in a variable (e.g.
% 'model_signal.mat') that you can load into the workspace everytime 
% after that

try 
    load('model_signal.mat')
catch
    figure(fig_handle);
    fprintf('\n Select anchor points of an ECG heartbeat (using mouse left click)');
    [xd,yd] = ginput % choose an undetermined number of points
    xd = round(xd*ECG.fs) % round towards integer sample indices
    model_signal.x = xd;
    %use the selected 'anchor points' to construct a linear model
    
end

piecwise_index = model_signal.x;
% Plot of model
figure
plot((piecwise_index)/fs,data(piecwise_index),':');
xlabel('Time (s)')
ylabel('Voltage (mV)')
axis('tight')
title('ECG Ideal signal model');
%
% Interpolate model to contain same time points as segmented signal

y_model=[];
for i=1:length(piecwise_index)-1
    %a=[i; i+1]
    %vq3 = interp1([model_signal.x(i),model_signal.x(i+1)],[data(piecwise_index(i)),data(piecwise_index(i+1))],model_signal.x(i):model_signal.x(i+1))
    if i==1 
        vq3 = interp1([piecwise_index(i),piecwise_index(i+1)],[data(piecwise_index(i)),data(piecwise_index(i+1))],piecwise_index(i):piecwise_index(i+1));
    else
        vq3 = interp1([piecwise_index(i),piecwise_index(i+1)],[data(piecwise_index(i)),data(piecwise_index(i+1))],piecwise_index(i)+1:piecwise_index(i+1));
    end
    y_model=[y_model vq3];
end

% --- Select noise segments
% Here, 5 noise segments will be selected that will be use to estimate the
% PSD of the noise. Repeat the selection, if necessary. 
fprintf('Go to the figure of the ECG signal and\nzoom in on a relevant portion of the signal\nin which you can observe the noise.\n')
%figure(fig_handle) % this will bring you back to the figure with the ECG signal
%keyboard % this will pause the execution, giving you time to execute the tip below:
% TIP: First, zoom in on a part of the ECG signal, in which you 
% can view several segments before making the selection. When you are
% ready, press F5 or click 'Run' to continue

% once you have found 'good' segments, save them in a variable (e.g.
% 'segments_noise.mat') that you can load into the workspace everytime 
% after that

nsegs = 5;

try 
    load('segments_noise.mat')
catch
    ans = 1;
    while ans~=0 % repeat the selection if desired

        % select segments and store the end points in x1,y1
        figure(fig_handle);
        fprintf('\n Select %d pairs of start and end points of segments (using mouse left click) : ',nsegs);
        [xn,yn] = ginput(nsegs*2);

        xn = round(xn*ECG.fs); % round towards integer sample indices
            % (Try to understand why this is not necessary for y1?)      
        ans = input(' \n Do you want to select new segments? (input 1 for yes) [1/0] : '); % enter 0 if the current selection is satisfactory

    end
    save('segments_noise.mat','xn')
end


% specifiy number of zeros for zero_padding 
n_zeros= 200;

% --- Desired signal's PSD estimation, based on the model signal

% Zero padding the signal 
n_zeros_signal = n_zeros-length(y_model);
y_model= [y_model zeros(1,n_zeros_signal)];

% Power spectrum of model
[amp_model,psd_model_P,hz_model,psd_model] = FourierT(y_model,fs);

% Plot of model 
figure     
plot((1:length(y_model))/fs,y_model)
title('Interpolated model')
% --- Take noise segments from the signal (based on the selected points), and compute the noise PSD

piecwise_index_n = xn;


% initalize matrix for storing noise segments
noise_mat= zeros(nsegs,n_zeros);

% vector to store lengths of noise segments 
li=zeros(1,5);

j=1;
% Create matrix with each noise segment (DC component removed) 
for i=1:nsegs
    noise_interval = (data(piecwise_index_n(j):(piecwise_index_n(j+1))));
    noise_mean = mean(noise_interval);
    n_zeros_seg = n_zeros-length(noise_interval);
    li(1,i)= length(noise_interval); %store values needed for normalization of PSD
    noise_interval=[noise_interval -  noise_mean]; 
    noise_interval= [noise_interval zeros(1,n_zeros_seg)];
    noise_mat(i,:) = noise_interval;
    j = j+2;
end

%Caluclate PSD of each noise segment 
PSD_noise_matrix_n = zeros(nsegs,n_zeros);

for i=1:nsegs
    [amplitude_noise,PSD_noise_seg_P,hz_noise,PSD_noise_seg] = FourierT(noise_mat(i,:),fs);
    PSD_noise_matrix_n(i,:) = PSD_noise_seg*(n_zeros/li(i));
end

% Computing the average of all different segments 
PSD_noise = sum(PSD_noise_matrix_n,1)/nsegs;

% --- Computation of the Wiener filter

Wiener_filt = psd_model./(PSD_noise+psd_model);


% --- Apply the Wiener filter to the signal
% --- (use the provided function _wola_filt2019_ and read its documentation)

open('wola_filt2019')

filt_data = wola_filt2019(data',Wiener_filt');


figure
hold on 
plot(t,data)
plot(t,filt_data)
xlabel('Time (s)')
ylabel('Voltage (mV)')
legend('Original','filtered signal (Wiener filter)')
title('Original vs filtered signal (Wiener filter)')
     
% --- plot the required PSDs and spectra. 
% --- (use the code below to arrange the graph)

relspacing = 0.25;
nrows = 3; vspace = relspacing/(nrows+1); height = (1-(nrows+1)*vspace)/nrows;
relspacing = 0.15;
ncols = 3; hspace = relspacing/(ncols+1); width = (1-(ncols+1)*hspace)/ncols;
figure('name','Wiener filtering');
ax = [];
idx = 3; ax(1) = subplot('position',[ hspace+mod(idx-1,ncols)*(width+hspace) , vspace+(nrows-ceil(idx/ncols))*(height+vspace) , width , height ]);

    plot(t,data)
    xlabel('Time (s)')
    ylabel('Voltage (mV)')
    title('\color[rgb]{1,0,0} input ECG signal','fontsize',10)
    
idx = 6; subplot('position',[ hspace+mod(idx-1,ncols)*(width+hspace) , vspace+(nrows-ceil(idx/ncols))*(height+vspace) , width , height ])	
    plot(hz_model,10*log10(Wiener_filt(1:length(hz_model))))
    xlabel('Frequency (hz)')
    ylabel('Log frequency response (dB)')
    axis('tight')
    title('\color[rgb]{0,0.3,0.7} Wiener transfer function - W(f)','fontsize',10)
    
idx = 9; ax(2) = subplot('position',[ hspace+mod(idx-1,ncols)*(width+hspace) , vspace+(nrows-ceil(idx/ncols))*(height+vspace) , width , height ]);
    plot(t,filt_data)
    xlabel('Time (s)')
    ylabel('Voltage (mV)')
    title('\color[rgb]{0,0.8,0} output ECG signal','fontsize',10)
    
idx = 4; subplot('position',[ hspace+mod(idx-1,ncols)*(width+hspace) , vspace+(nrows-ceil(idx/ncols))*(height+vspace) , width , height ])
    plot(hz_model,10*log10(psd_model_P))
    grid on
    xlabel('Frequency (Hz)')
    ylabel('PSD (in dB)')
    axis tight;
    title('PSD of desired signal -  S_d(f)','fontsize',10)
    
idx = 5; subplot('position',[ hspace+mod(idx-1,ncols)*(width+hspace) , vspace+(nrows-ceil(idx/ncols))*(height+vspace) , width , height ])   	
    plot(hz_model,10*log10(PSD_noise(1:length(hz_model))))
    grid on
    xlabel('Frequency (Hz)')
    ylabel('PSD (in dB)')
    xlim([1.25 hz_model(end)])
 
    title('PSD of noise - S_n(f)','fontsize',10)
linkaxes(ax)    

    figure
    plot(hz_model,10*log10(PSD_noise(1:length(hz_model))))
    xlim([20 30])
    grid on
    xlabel('Frequency (Hz)')
    ylabel('PSD (in dB)')
    title('PSD of noise - S_n(f)','fontsize',10)

