function [] = topoPlot(A)
% TOPOPLOT Plot a topographic map.
%   TOPOPLOT(A) plots a topographic maps based on the matrix A of
%   contributions of each component (e.g., a mixing matrix).

sz = 50;
[T,crds] = topoTensor(A,sz);

for i = 1:size(T,3)

    subplot(4,3,i);
    
    im = T(:,:,i);
    imagesc(abs(im)); hold on; colormap jet;
    
    d = size(im,1);
    t = 0:pi/128:(2*pi);
    cx = (d-1)/2*cos(t)+(d+1)/2;
    cy = (d-1)/2*sin(t)+(d+1)/2;
    plot(cx,cy,'k--','lineWidth',2);
    scatter(crds(:,1),crds(:,2),50,'k','filled');
    axis off
    
    title(['Component ' num2str(i)]);

end
    
end

function [topoTens, crds] = topoTensor(EEGmat, imsize)
%[topoTens, crds] = topoTensor(EEGmat, imsize)
%
%Calculates a topograpic tensor from a channel-time EEG matrix. The
%spatial dimension sizes are defined as imsize, the temporal dimension size
%equals the number of columns in the EEG matrix
%
% EEGmat: a channel x time matrix of EEG data
% imsize: the side length of the square topographic map for one time instant
% system: the electrode system ('1020', 'debener14', 'debener12');
%           1020 is the standard 10-20 system
%           debener14 has the fourteen channels of Debener's input cap
%           debener12 is the same, but Tp9 and Tp10 are eliminated through rereferencing
% topoTens: the constructed topographic tensor (spatial1 x spatial2 x time)
% crds: coordinates of the electrodes in image coordinates

% Billiet 2013-2016

system = {'O2','O1','P4','Pz','P3', 'C4' ,'Cz','C3','Fz','F4','F3','Fpz'};

load EEGcoordsystem.mat;

[crds, sigmas] = calcCoordinatesAndSigmas(imsize, system, coordlabels, coordsmat);

topoTens = zeros(imsize,imsize,size(EEGmat,2));
for imInd = 1:size(EEGmat,2)
    topoTens(:,:,imInd) = calcSingleMap(imsize, crds, sigmas, EEGmat(:,imInd));
end

end

function im = calcSingleMap(imsize, crds, sigmas, EEGsample)
%calculates a single topographic map from given size, coordinates, gaussian
%widths and the EEG sample values for each electrode. The map is based on
%combinning Gaussian distributions on the electrode positions.

% Billiet 2013-2015

    all = zeros(imsize,imsize,length(EEGsample));
    [X Y] = meshgrid(1:imsize, 1:imsize);
    wMat = zeros(length(EEGsample),length(EEGsample));
    for i = 1:length(EEGsample)
       imComp = mvnpdf([X(:) Y(:)],crds(i,1:2),sigmas(i)*eye(2));
       wMat(:,i) = EEGsample(i)*mvnpdf(crds(:,1:2),crds(i,1:2),sigmas(i)*eye(2));
       all(:,:,i) = reshape(imComp,imsize,imsize)*EEGsample(i);
    end
    weights = wMat\EEGsample;
    for j = 1:length(EEGsample)
       all(:,:,j) = weights(j)*all(:,:,j); 
    end
    im = (sum(all,3));
end

function [crds, sigmas] = calcCoordinatesAndSigmas(imsize, system, coordlabels, coordsmat)
% Calculates the electrodes' coordinates based on the image size and the
% coordinate system, as well as the widths of the gaussian distributions
% that will represent their activity

% Billiet, 2013-2015

l = numel(system);
crds = zeros(l, 3);
for ic = 1:l
    ind = find(strcmpi(system{ic},coordlabels));
    crds(ic,:) = coordsmat(ind,:);
end

crds = crds*(imsize-1)+1; %remap from 0:1 to 1:imsize
temp = crds;
crds(:,1) = imsize+1 - temp(:,2);
crds(:,2) = imsize+1 - temp(:,1);

d = dist(crds(:,1:2)');
d(d==0) = Inf;
sigmas = (min(d)/1.75).^2;
end