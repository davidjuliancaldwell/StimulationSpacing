%% script to analyze stimulation pulses for subject 2fd831

close all;clear all;clc

sid = '2a8d14';
SUB_DIR = fullfile(myGetenv('subject_dir'));
montageFilepath = strcat(SUB_DIR,'\',sid,'\',sid,'_Montage.mat');
load(montageFilepath);
load('2a8d14_mappedElectrodes.mat')
filePath = 'C:\Users\djcald.CSENETID\GoogleDrive\GRIDLabDavidShared\2a8d14\StimulationSpacingChunked\stim_widePulse_9_33_15_21';
load(filePath)


%% plot individual trials for each condition on a different graph

labels = unique(abs(Sing(Sing~=0)));
uniqueLabels = unique(abs(Sing(Sing~=0)));

% intialize counter for plotting
k = 1;

% make vector of stim channels
stimChans = stim_chans;

% determine number of subplot
numChans = size(dataEpoched,2);
subPlots = numSubplots(numChans);
p = subPlots(1);
q = subPlots(2);
scaling = false;

realStim = [9 33];

% plot each condition separately e.g. 1000 uA, 2000 uA, and so on

for i=uniqueLabels
    figure;
    dataInterest = dataEpoched(:,:,:);
    for j = 1:numChans
        subplot(p,q,j);
        plot(t,squeeze(dataInterest(:,j,:)));
        xlim([min(t) max(t)]);
        
        % change y axis scaling if necessary
        if strcmp(scaling,'y')
            ylim([minVal maxVal]);
        end
        
        % put a box around the stimulation channels of interest if need be
        if ismember(j,stimChans)
            ax = gca;
            ax.Box = 'on';
            ax.XColor = 'red';
            ax.YColor = 'red';
            ax.LineWidth = 2;
            title(num2str(j),'color','red');
        elseif ismember(j,realStim)
            ax = gca;
            ax.Box = 'on';
            ax.XColor = 'blue';
            ax.YColor = 'blue';
            ax.LineWidth = 2;
            title(num2str(j),'color','blue');
        else
            title(num2str(j));
            
        end
        vline(0);
        
    end
    
    % label axis
    xlabel('time in ms');
    ylabel('voltage in V');
    subtitle(['Individual traces - Current set to ',num2str(uniqueLabels(k)),' \muA']);
    
    
    % get cell of raw values, can use this to analyze later
    dataRaw{k} = dataInterest;
    
    % get averages to plot against each for later
    % cell function, can use this to analyze later
    dataAvgs{k} = mean(dataInterest,3);
    
    
    
    %increment counter
    k = k + 1;
    
    
end

%% Plot Cortex 
badChans = [9 15 21 29 30 33];
% colorbrewer colormap
% flip it so red is increase, blue is down
CT = cbrewer('seq','PuRd',128);

stimChans = [9 33];
figure

hold on
gChansMat = logical(ones(size(locs,1),1));
gChansMat(badChans) = 0;
gChansLabel = [1:size(locs,1)];
gChansLabel = gChansLabel(gChansMat);
locsG = locs(gChansMat,:);
weights = max(mean(dataEpoched(:,gChansMat,:),3));

PlotDotsDirect(sid,locsG,weights,'b',[0 abs(max(weights))],25,CT,gChansLabel,true,false,0.25);

weightsNew = ones(size(stimChans));
locsC = locs(stimChans,:);
PlotDotsDirect(sid,locsC,weightsNew,'b',[0 abs(max(weights))],25,[0.1 0.5 1; 0.1 0.5 1;0.1 0.5 1],stimChans,true,true,0.25);
 colormap(CT)
cbar = colorbar;
ylabel(cbar,'max average stimulation pulse magnitude (V)')
title('Peak Stimulation Voltage')
% stim chans

%%

% CT = flipud(CT);
% CT = cbrewer('seq','PuRd',128);
% stimChans = [9 33];
% 
% figure
% 
% hold on
% weights = max(mean(dataEpoched(:,:,:),3));
% weights(stimChans) = 0;
% 
% 
% PlotDotsDirect(sid,locs,weights,'b',[0 abs(max(weights))],15,CT,1:128,true,false,0.25);

