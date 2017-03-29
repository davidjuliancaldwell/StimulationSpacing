%% script to analyze stimulation pulses for subject 2fd831

sid = '2fd831';
SUB_DIR = fullfile(myGetenv('subject_dir'));
montageFilepath = strcat(SUB_DIR,'\',sid,'\',sid,'_Montage.mat');
load(montageFilepath);
load('2fd831_mappedElectrodes.mat')
filePath = 'C:\Users\djcald.CSENETID\GoogleDrive\GRIDLabDavidShared\2fd831\StimulationSpacingChunked\stim_widePulse_121_55';
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

weights = max(mean(dataEpoched,3));

% colorbrewer colormap
CT = cbrewer('div','RdBu',128);
% flip it so red is increase, blue is down

CT = flipud(CT);
CT = cbrewer('seq','PuRd',128);
figure
hold on

PlotDotsDirect(sid,locs,weights,'b',[0 abs(max(weights))],15,CT,1:128,false,false,0.25);

colormap(CT)
colorbar;