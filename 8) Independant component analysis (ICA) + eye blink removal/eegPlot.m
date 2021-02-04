function eegPlot(S,meas,sens,linestyle)
% EEGPLOT Plot a multichannel EEG signal.
%   EEGPLOT(S,meas,sens,linestyle) plots a multichannel EEG signal S, 
%   containing m rows (channels) and n columns (time samples). Meas
%   represents the montage (a cell array of electrode names), sens is the
%   y-axis scale (optional) and there is a linestyle argument (optional).

if nargin<2, mode = 'raw'; elseif isempty(meas), mode = 'raw'; else mode = 'montage'; end

if nargin<4, linestyle='b'; end

fs = 128; % sample rate

[r c] = size(S);

if strcmp(mode,'montage')
  
    chanlabels = meas;
    nofchan = r;
  
elseif strcmp(mode,'raw')

  chanlabels=[];
  nofchan = r;
  for i=1:nofchan
    chanlabels{i} = num2str(i);
  end

end
%figure
chanlabels = chanlabels';
S=flipud(S);

[r c]=size(S);
if nargin<3
 m = max(max(abs(S)));
else 
 m = sens*2048/2000;
end
% m = 200;

plot([0,c/fs],[1 1],'k:')
% hold on
set(gca,'nextplot','add');
for i=1:r
  %m = max(abs(S(i,:)));
  if (m~=0)
    S(i,:)=S(i,:)/m + i;
  else
    S(i,:)=S(i,:) + i;  
  end
  plot([0,c/fs],[i i],'k:');
end
plot((1:c)/fs,S',linestyle,'linewidth',0.25,'color','k');
set(gca,'xlim',[0 c/fs]);
set(gca,'ylim',[0 r+1]);
set(gca,'ytick',1:nofchan);
set(gca,'yticklabel',flipud(chanlabels));
xlabel('Time (sec)')
if nargin<3
  sensstring=num2str(round(m*2000/2048));
else
  sensstring=num2str(sens);
end

text(0,r+0.5,[sensstring ' uV'])
hold off
