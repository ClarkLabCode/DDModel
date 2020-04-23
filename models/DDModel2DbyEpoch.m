function a = DDModel2DbyEpoch(stimPath,varargin)

%% Created by Ryosuke Tanaka (cleaned up for publication 4/23/20)
% This function performs a numerical simulation of the
% displacement-detector model of Drosophila LC11 reported in Tanaka & Clark
% (2020) Curr. Biol. 
% The funciton takes a path to folders that include visual stimuli (XYT
% plots) and their metadata, and returns a structure with simulated
% responses of the model.

% The overall model architecture is as follows:
% 1. Half-wave rectification (modeling ON-OFF responses)
% 2. Center-surround antagonism (modeling size tuning)
% 3. Feedforward inhibition (adaptation)
% 4. Spatiotemporal pooling (LC11 voltage)

%% Default parameters
% Ones concerning the interpretation of XYT plot

% XYT plot files are essentially a large csv, whose first 5 columns are
% metadata and the rest contains actual visual stimuli. The first 2 columns
% are timestamps unrelated to the purpose here. Third column indicates which 
% stimulus was presented (1 corresponds to interleaves). The function
% assuems that interleave duration does not change.
% Forth and fifth columns encode the spatial dimension of the stimulus.
% This is used to recover flattened stimulus.
% The last column is always empty because of the way those csv is
% generated.

% The stimuli were originally generated at the resolution of 1 deg for
% physiological experiment, and then downsampled using interpolation.
% We did not include the stimulus generation code since these are far
% too complicated for the purpose here (because it was designed for
% on-line stimulus presentation for experiments), but it should be straight-
% forward for anyone interested to generate a 180Hz/5deg resolution stimulus
% as a 3D matrix in matlab, flatten it, append necessary meta-data to feed
% into this function.

% We assume XYT plot to be already in the ommatidial resolution of 5
% degrees. Temporal resolution is 180Hz.
frameRate = 180; % (Hz)
XYTbitDepth = 4;   % bit depth of stimulus (this is set to 4 to keep file size small)
meanGray = ceil((2^XYTbitDepth - 1)/2);
trimXYT = [22,22]; % trim the stimulus (defaulted at 22 ommatidia = 110 deg)

% Model parameters
% initial temporal high-pass filter
tau_pre = 200; % (ms)

% center-surround antagonism parameters
kerSizes = [5,15];    % Sigma of Gaussian for the excitatory and inhibitory
                      % lobes of the upstream spatial filtering, in degree.
                      % Correponds to sigma1 and sigma2 in the paper.
                      
kerWeights = [1,3.5]; % Weights of the excitatory and inhibitory lobes.
                      % Corresponds to w1 and w2 in the paper.
                      
filterOrder = 2;      % for non-Gaussian spatial kernel 

% adaptation parameters
tau_q = 300;    % tau of adaptation state (ms, tau_adapt in the paper)
alpha = 1;      % gain of stimulus into response
gamma = 1000;   % gain of adaptation state onto response

% spatial pooling
sigma = 10;  % (deg) sigma_out in the paper
tau   = 300; % (ms)  tau_out in the paper

% how much duration to simulate (around stimulus onset)
snipDuration = 10; % (s) entire duration of each epoch
snipShift = -2;    % (s) pre-stimulus period

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

% Read out the spatial dimension of the stimulus
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
nStep = 6;

%% Do simulation
% prepare output structure
a = struct;
a.timeX = timeX;
a.snipShift = snipShift;
figure;
for ee = 1:numEpochs
    
    % Cut stimulus snippet out of the entire xytplot, fill the rest with gray
    epochStart = sum(epochFrames(1:ee-1)) + ILframeCount*ee;
    epochEnd = sum(epochFrames(1:ee)) + ILframeCount*ee;
    stimSnip = xytp(epochStart+1:epochEnd,:);
    S = meanGray*ones(length(timeX),size(xytp,2));
    S((1:epochFrames(ee)) - snipShift*frameRate,:) = stimSnip;
    
    % convert luminance to contrast
    S = (S-meanGray)/meanGray;

    % convert the flattened stimulus back into space 2D
    S = permute(reshape(S', [dims size(S,1)]),[2,1,3]);

    % trim the stimulus in space
    if any(trimXYT)
        S = S((1:trimXYT(1))+round((size(S,1)-trimXYT(1))/2),...
              (1:trimXYT(2))+round((size(S,2)-trimXYT(2))/2),:);
    end
    sizeX = size(S,2);
    sizeY = size(S,1);
    
    % show stimulus time trace (on a horizontal slice in the middle)
    subplot(numEpochs,nStep,nStep*(ee-1)+1);
    plot(permute(S(round(sizeX/2),:,:),[3,2,1])); title('Stimulus');
    drawnow;

    %% 0. High-pass filterling
    disp('High pass filtering...');
    kerTpre = (tau_pre-timeX).* exp(-timeX/tau_pre)/tau_pre/tau_pre;
    kerTpre = kerTpre / sqrt(sum(kerTpre.^2)/frameRate); % unit l2 norm  
    kerTpre(timeX>tau_pre*10) = [];
    RHP = convn(S,permute(kerTpre,[3,1,2]),'full')/frameRate;
    RHP(:,:,length(timeX)+1:end) = [];
    subplot(numEpochs,nStep,nStep*(ee-1)+2); plot(permute(RHP(round(sizeX/2),:,:),[3,2,1])); title('High-pass filtering');
    drawnow;
    
    %% 1. Fullwave rectification  %%
    Rrect  = RHP .*(RHP>0) - RHP.*(RHP<0); 

    %% 2. Center surround %%
    disp('Spatial filtering...');
    RCS  = DoGrect(Rrect,kerSizes,kerWeights,filterOrder); % subroutine for this (at the end)
    subplot(numEpochs,nStep,nStep*(ee-1)+3); 
    plot(permute(RCS(round(sizeX/2),:,:),[3,2,1])); 
    title('After center-surround');

    %% 3. Adaptation %%
    % Divide RCS with its low-pass version
    kerTadapt = timeX.* exp(-timeX/tau_q)/tau_q/tau_q;
    kerTadapt = kerTadapt / sqrt(sum(kerTadapt.^2)/frameRate); % unit l2 norm
    
    Q = convn(RCS,permute(kerTadapt,[3,1,2]),'full')/frameRate;
    Q = Q(:,:,1:size(RCS,3)); % truncate
    Radapt = alpha*RCS./(1+gamma*Q);
    subplot(numEpochs,nStep,nStep*(ee-1)+4); plot(permute(Q(round(sizeX/2),:,:),[3,2,1])); title('Q');
    subplot(numEpochs,nStep,nStep*(ee-1)+5); plot(permute(Radapt(round(sizeX/2),:,:),[3,2,1])); title('R');

    %% 4. Spatial pooling %%
    % separable spatio temporal low-pass filter 
    kerXSize = ceil(sigma/5)*5*3; % spatial kernal size = 3 sigma 
    degs = -kerXSize:5:kerXSize;
    kerX = exp(-degs.^2/2/sigma^2)/sqrt(2*pi*sigma^2);
    kerT = timeX .* exp(-timeX/tau)/tau/tau;
    kerT = kerT / sqrt(sum(kerT.^2)/frameRate); % unit l2 norm
    kerT(timeX>tau*10) = [];

    Rout = convn([Radapt, Radapt, Radapt],kerX,'same')*5;
    Rout = convn([Rout; Rout; Rout],kerX','same')*5;
    Rout = Rout((1:sizeY)+sizeY,(1:sizeX)+sizeX,:)/frameRate;

    Rout = convn(Rout,permute(kerT,[3,1,2]),'full');
    Rout(:,:,length(timeX)+1:end) = [];
    
    subplot(numEpochs,nStep,nStep*(ee-1)+6); 
    plot(permute(Rout(round(sizeX/2),:,:),[3,2,1])); 
    title('Output');

    % record results
    a(ee).S = S;
    a(ee).RHP = RHP;
    a(ee).Rrect = Rrect;
    a(ee).RCS = RCS;
    a(ee).Q = Q;
    a(ee).Radapt = Radapt;
    a(ee).Rout = Rout;
end

%% Visualization

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


