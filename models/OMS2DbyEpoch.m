function a = OMS2DbyEpoch(stimPath,varargin)

%% Created by Ryosuke Tanaka (cleaned up for publication 4/23/20)
% This function performs a numerical simulation of the inprementation of
% OMS model of Drosophila LC11 reported in Tanaka & Clark
% (2020) Curr. Biol. The model architecutre is based on Olveczky et al.
% (2003) Nature.
% The funciton takes a path to folders that include visual stimuli (XYT
% plots) and their metadata, and returns a structure with simulated
% responses of the model.

% The overall architecture is as follows:
% 1. Front end temporal high-pass filtering
% 2. rectification
% 3. Center surround antagonism

% Chop up the stimulus matrix into epochs so that accumulating long term
% error doesn't affect the model

%% Default parameters
% Ones concerning the interpretation of XYT plot
% See details for DDModel2DbyEpoch
XYTbitDepth = 4;
meanGray = ceil((2^XYTbitDepth - 1)/2);
trimXYT = [22,22];
frameRate = 180;
% initial temporal high-pass filter
tau_pre = 200;

% center-surround parameters
kerSizes = [7.5,20];
kerWeights = [1,3.5];
filterOrder = 2;

% calcium filter (low pass)
tau = 300;

% how much duration to simulate (around stim onset)
snipDuration = 10;
snipShift = -2;


% for overwriting any of default parameters
for ii = 1:2:length(varargin)
    eval([varargin{ii} '= varargin{' num2str(ii+1) '};']);
end
%% Data preparation
% open GUI to get path for data if it's not provided
thisdir = pwd; % assuming this is run in the DDModel folder
if nargin<1
    stimPath = [];
end

if isempty(stimPath)
    DDmodelstrind = strfind(thisdir,'DDModel');
    if ~isempty(DDmodelstrind)
        stimDir = [thisdir(1:DDmodelstrind+6),'/stimuli',];
    else
        stimDir = [];
    end
    stimPath = uigetdir(stimDir);
end
xytp = csvread([stimPath,'/xtPlot.xtp']); % read stimulus
load([stimPath,'/epochNames.mat']);

dims = [xytp(1,4),xytp(1,5)];
% Read out stimulus type timetrace
epochTrace = xytp(:,3);
epochOnset  = [0;diff(epochTrace)>0];  % onset of epochs
epochOffset = [0;diff(epochTrace)<0];% offset of epochs

% discard non stimulus columns
xytp = xytp(:,6:end-1);

% Identify the length of interleave (i.e. length of the first epoch)
ILframeCount = find(epochOnset==1,1)-1;

% find out how many epochs there are
numEpochs = max(epochTrace)-1;
% find out the length of each epoch
epochFrames = find(epochOffset==1)-find(epochOnset==1);

% frameRate = 60*params(1).framesPerUp;
% prepare a vector with time steps
timeX = linspace(1,snipDuration*1000,snipDuration*frameRate);

% number of steps to visualize
nStep = 4;


%% Do simulation
a = struct;
a.timeX = timeX;
a.snipShift = snipShift;
figure;
for ee = 1:numEpochs
    
    % Cut stimulus snippet out of the entire xytplot, fill with gray
    epochStart = sum(epochFrames(1:ee-1)) + ILframeCount*ee;
    epochEnd = sum(epochFrames(1:ee)) + ILframeCount*ee;
    stimSnip = xytp(epochStart+1:epochEnd,:);
    S = meanGray*ones(length(timeX),size(xytp,2));
    S((1:epochFrames(ee)) - snipShift*frameRate,:) = stimSnip;
    
    % convert to contrast
    S = (S-meanGray)/meanGray;

    % convert the image into ommatidial resolution (5deg)
    S = permute(reshape(S', [dims size(S,1)]),[2,1,3]);

    if any(trimXYT)
        S = S((1:trimXYT(1))+round((size(S,1)-trimXYT(1))/2),...
              (1:trimXYT(2))+round((size(S,2)-trimXYT(2))/2),:);
    end
    sizeX = size(S,2);
    sizeY = size(S,1);
    sizeT = size(S,3);
    % show stimulus
    subplot(numEpochs,nStep,nStep*(ee-1)+1);
    plot(permute(S(round(sizeX/2),:,:),[3,2,1]));
    title('Stimulus');
    drawnow;

    %% 1. High-pass filterling
    kerTpre = (tau_pre-timeX).* exp(-timeX/tau_pre)/tau_pre/tau_pre;
    kerTpre = kerTpre / sqrt(sum(kerTpre.^2)/frameRate); % unit l2 norm  
    kerTpre(timeX>tau_pre*10) = [];
    RHP = convn(S,permute(kerTpre,[3,1,2]),'full')/frameRate;
    RHP(:,:,length(timeX)+1:end) = [];
    subplot(numEpochs,nStep,nStep*(ee-1)+2);
    plot(permute(RHP(round(sizeX/2),:,:),[3,2,1]));
    title('High-pass filtering');
    drawnow;
    
    %% 2. Halfwave rectification  & size tuning %%
    RON  = DoGrect(RHP.*(RHP>0) ,kerSizes,kerWeights,filterOrder);
    
    subplot(numEpochs,nStep,nStep*(ee-1)+3);
    plot(permute(RON(round(sizeX/2),:,:),[3,2,1]));
    title('Surround supression (ON)');
    
    %% 3. Rear-end filtering %%
    kerTout = timeX.* exp(-timeX/tau)/tau/tau;
    kerTout = kerTout / sqrt(sum(kerTout.^2)/frameRate); % unit l2 norm
    Rout = convn(RON,permute(kerTout,[3,1,2]),'full')/frameRate;
    Rout = Rout(:,:,1:size(RON,3));
    subplot(numEpochs,nStep,nStep*(ee-1)+4);
    plot(permute(Rout(round(sizeX/2),:,:),[3,2,1]));
    
    % record results
    a(ee).S = S;
    a(ee).RHP = RHP;
    a(ee).RON = RON;
    a(ee).Rout = Rout;
end


%% Summary Visualization

% time trace
figure; subplot(1,2,1);
meanResponses = nan(numEpochs,1);
for ee = 1:numEpochs
    plot(permute(a(ee).Rout(round(sizeY/2),round(sizeX/2),:),[3,2,1])); hold on
    meanResponses(ee) = mean(a(ee).Rout(round(sizeY/2),round(sizeX/2),:),3);
end
legend(epochNames);

% integrated responses
subplot(1,2,2); bar(meanResponses);
xticks(1:numEpochs);
xticklabels(epochNames);

end

function out = DoGrect(M,sigs,ws,filterOrder)
% DoG filtering followed by rectification
% for within-polarity antagonism
if nargin<4
    filterOrder = 2;
end
if nargin<3
    ws = [1,1];
end
kerSize = ceil(max(sigs)/5)*5*3;
[degsX,degsY] = meshgrid(-kerSize:5:kerSize);
gauss1 = exp(-(degsX.^filterOrder+degsY.^filterOrder)/2/sigs(1)^filterOrder)/(2*pi*sigs(1)^filterOrder);
gauss2 = exp(-(degsX.^filterOrder+degsY.^filterOrder)/2/sigs(2)^filterOrder)/(2*pi*sigs(2)^filterOrder);
ker    = gauss1*ws(1) - gauss2*ws(2); 
ker    = ker ./ sqrt(sum(ker(:).^2)*5*5); % unit l2 norm


% padding
M2 = [M(end:-1:1,end:-1:1,:),M(:,end:-1:1,:),M(end:-1:1,end:-1:1,:);...
      M(:,end:-1:1,:),M,M(:,end:-1:1,:);...
      M(end:-1:1,end:-1:1,:),M(:,end:-1:1,:),M(end:-1:1,end:-1:1,:)];
out = convn(M2,ker,'same')*5*5;
out = out(size(M,1)+(1:size(M,1)),size(M,2)+(1:size(M,2)),:);
out = out.*(out>0);
end


