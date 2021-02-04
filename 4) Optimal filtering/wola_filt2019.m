function [ y ] = wola_filt2019( x , H )
%WOLA_FILT implements Weighted Overlap-Add Filtering, for sliding window
%filtering of long signals.
%
%   INPUTS
%   x = matrix with input signals to be filtered (samples x variables)
%   H = filter in frequency domain (frequencies x 1)
%
%   OUTPUTS
%   y = matrix with filtered output signals (samples x variables)

%% Initialize
% -- Determine sizes
sig_length = size(x,1);
nsigs = size(x,2);
L = length(H); % FFT length ( = filter length)
    assert(mod(L,2)==0)
M = 1*L; % window length
    assert(M<=L)
sh = round(M/2); % window shift
nw = floor((sig_length-L)/sh)+1; % number of windows

% -- Create window function
win = 'hann';
if ~exist('win','var')
    w = window(@hann,M);
else
    assert(ischar(win))
    w = eval(strcat('window(@',win,',M)'));
end
w = sqrt(w(:));

% -- Apply zero-padding if the FFT length and window length differ
if L > M
    w = [ w ; zeros(L-M,1) ]; % pad with zeros
end

%% Analysis + filtering + synthesis
y = zeros((nw-1)*sh+L,nsigs);
x_window = zeros(L,nsigs); x_weighted = zeros(L,nsigs);
y_window = zeros(L,nsigs); y_weighted = zeros(L,nsigs);
for i = 1:nw
    start = (i-1)*sh+1;
    stop = start+L-1;
    % analysis
    x_window = x(start:stop,:);
    x_weighted = bsxfun(@times,w,x_window);
    X = fft(x_weighted);
    % filtering
    Y = bsxfun(@times,X,H);
    % synthesis
    y_window = ifft(Y);
    y_weighted = bsxfun(@times,w,y_window);
    y(start:stop,:) = y(start:stop,:) + y_weighted;
end   

%% Pad the output if nececssary, to match the input data length
if size(y,1) < sig_length
    y = [ y ; zeros(sig_length-size(y,1),nsigs) ];
end

end

