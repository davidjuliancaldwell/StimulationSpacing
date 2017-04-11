%% 6-23-2016 - David Caldwell - script to look at stim spacing
% starting with subject 3f2113

%% initialize output and meta dir
% clear workspace
close all; clear all; clc

% set input output working directories
Z_ConstantsStimSpacing;
%CODE_DIR = fullfile(myGetenv('gridlab_dir'));
%scripts_path =  strcat(CODE_DIR,'\Experiment\BetaTriggeredStim\scripts');

% add path for scripts to work with data tanks
addpath('../scripts')

% subject directory, change as needed
SUB_DIR = fullfile(myGetenv('subject_dir'));

%% load in subject

% this is from my z_constants

sid = SIDS{4};

% load in tank
if (strcmp(sid, '3f2113'))
    tank = TTank;
    tankSelect = strcat(SUB_DIR,'\',sid,'\data\data\d6\stimMoving\stimMoving');
    tank.openTank(tankSelect);
    
    % select the block
    tank.selectBlock('stimSpacing-22');
    %  mark stim channels if desired
    stim_chans = input('Input the stim channels as an array e.g. [22 30]');
    % stims = [29 28];
    
    % load in the data, stim info, sampling rates
    tic;
    [data, data_info] = tank.readWaveEvent('Wave');
    [stim, stim_info] = tank.readWaveEvent('Stim');
    toc;
    
    % get sampling rates
    fs_data = data_info.SamplingRateHz;
    fs_stim = stim_info.SamplingRateHz;
    
    % get current delivery
    [Stm0, Stm0_info] = tank.readWaveEvent('Stm0');
    [Sing, Sing_info] = tank.readWaveEvent('Sing');
    % elseif (strcmp(sid,'20f8a3'))
    %     tank = TTank;
    %     tankSelect = strcat(SUB_DIR,'\',sid,'\data\d7\SensoryScreen\sensoryscreen');
    %     tank.openTank(tankSelect);
    %
    %     % select the block
    %     tank.selectBlock('sensoryScreen-8');
    %     %  mark stim channels if desired
    %     stim_chans = input('Input the stim channels as an array e.g. [22 30]');
    %     % stims = [29 28];
    %
    %     % load in the data, stim info, sampling rates
    %     tic;
    %     [data1, data_info] = tank.readWaveEvent('ECO1');
    %     [data2, data_info] = tank.readWaveEvent('ECO2');
    %
    %     [data3, data_info] = tank.readWaveEvent('ECO3');
    %
    %     [data4, data_info] = tank.readWaveEvent('ECO4');
    %
    %
    %     [stim, stim_info] = tank.readWaveEvent('Stim');
    %     toc;
    %
    %     % get sampling rates
    %     fs_data = data_info.SamplingRateHz;
    %     fs_stim = stim_info.SamplingRateHz;
    %
    %     % get current delivery
    %     [Stm0, Stm0_info] = tank.readWaveEvent('Stm0');
    %     [Sing, Sing_info] = tank.readWaveEvent('Sing');
    %
    %     data = [data1 data2 data3 data4];
    %     clear data1 data2 data3 data4
    %
    %     OUTPUT_DIR = 'C:\Users\djcald.CSENETID\Google Drive\GRIDLabDavidShared\20f8a3\StimulationSpacingChunked';
    
    %%%%%%%%% below is for wide!!!
elseif (strcmp(sid,'20f8a3'))
    filePath = 'C:\Users\djcald.CSENETID\GoogleDrive\GRIDLabDavidShared\20f8a3\StimulationSpacing\StimSpacing_28_4_Wide';
    
    prompt = {'how many channels did we record from? e.g 48 ', 'what were the stimulation channels? e.g 28 29 ', 'how long before each stimulation do you want to look? in ms e.g. 1', 'how long after each stimulation do you want to look? in ms e.g 5'};
    dlg_title = 'StimChans';
    num_lines = 1;
    defaultans = {'128','28 4','1','10'};
    answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
    numChans = str2num(answer{1});
    stim_chans = str2num(answer{2});
    preTime = str2num(answer{3});
    postTime = str2num(answer{4});
    
    stim_1 = stim_chans(1);
    stim_2 = stim_chans(2);
    
    %totalFile = strcat(filePath,num2str(stim_1),'_',num2str(stim_2),'inputFile);
    load(filePath);
    
    
    %%
    % ui box for input for stimulation channels
    
    % get sampling rates
    fs_data = Wave.info.SamplingRateHz;
    fs_stim = Stim.info.SamplingRateHz;
    
    % stim data
    stim = Stim.data;
    stim_info = Stim.info;
    % current data
    Sing_info = Sing.info;
    
    Sing = Sing.data;
    
    % recording data
    data = Wave.data;
    data_info = Wave.info;
    
    OUTPUT_DIR = 'C:\Users\djcald.CSENETID\GoogleDrive\GRIDLabDavidShared\20f8a3\StimulationSpacingChunked';
    
    
    
elseif (strcmp(sid,'2fd831'))
    %%
    filePath = strcat(SUB_DIR,'\',sid,'\data\d7\ConvertedMatlabFiles\stimSpacing\stimSpacing-');
    inputFile = input('which block ? \n','s');
    totalFile = strcat(filePath,inputFile);
    load(totalFile);
    
    Wave.info = ECO1.info;
    
    % ui box for input for stimulation channels
    prompt = {'how many channels did we record from? e.g 48 ', 'what were the stimulation channels? e.g 28 29 ', 'how long before each stimulation do you want to look? in ms e.g. 1', 'how long after each stimulation do you want to look? in ms e.g 5'};
    dlg_title = 'StimChans';
    num_lines = 1;
    defaultans = {'128','121 122','1','10'};
    answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
    numChans = str2num(answer{1});
    chans = str2num(answer{2});
    preTime = str2num(answer{3});
    postTime = str2num(answer{4});
    
    stim_1 = chans(1);
    stim_2 = chans(2);
    
    stim_chans = chans;
    
    % get sampling rates
    fs_data = Wave.info.SamplingRateHz;
    fs_stim = Stim.info.SamplingRateHz;
    
    % stim data
    stim = Stim.data;
    stim_info = Stim.info;
    % current data
    sing = Sing.data;
    Sing_info = Sing.info;
    Sing_info = Sing.info;
    
    Sing = Sing.data;
    
    % recording data
    data = [ECO1.data ECO2.data ECO3.data ECO4.data];
    data_info = Wave.info;
    
    
    OUTPUT_DIR = 'C:\Users\djcald.CSENETID\GoogleDrive\GRIDLabDavidShared\2fd831\StimulationSpacingChunked';
elseif (strcmp(sid,'2a8d14'))
    %%

   % ui box for input for stimulation channels
    prompt = {'how many channels did we record from? e.g 48 ', 'what were the stimulation channels? e.g 28 29 ', 'how long before each stimulation do you want to look? in ms e.g. 1', 'how long after each stimulation do you want to look? in ms e.g 5'};
    dlg_title = 'StimChans';
    num_lines = 1;
    defaultans = {'40','15 36','1','10'};
    answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
    numChans = str2num(answer{1});
    chans = str2num(answer{2});
    preTime = str2num(answer{3});
    postTime = str2num(answer{4});
    
    stim_1 = chans(1);
    stim_2 = chans(2);
    
    filePath = 'C:\Users\djcald.CSENETID\Data\ConvertedTDTfiles\2a8d14\stimSpacing_';
    
    
    totalFile = strcat(filePath,num2str(stim_1),'_',num2str(stim_2),'_wide');
    load(totalFile);
    
    
    
    stim_chans = chans;
    
    % get sampling rates
    fs_data = Wave.info.SamplingRateHz;
    fs_stim = Stim.info.SamplingRateHz;
    
    % stim data
    stim = Stim.data;
    stim_info = Stim.info;
    % current data
    sing = Sing.data;
    Sing_info = Sing.info;
    Sing_info = Sing.info;
    
    Sing = Sing.data;
    
    % recording data
    data = Wave.data;
    data_info = Wave.info;
    
    dataTemp = data(:,1:numChans);
    clear data;
    data = dataTemp;
    clear dataTemp;
    OUTPUT_DIR = 'C:\Users\djcald.CSENETID\GoogleDrive\GRIDLabDavidShared\2a8d14\StimulationSpacingChunked';

end
%%
plotIt = 'y';

%% plot stim
%
% figure
% hold on
% for i = 1:size(stim,2)
% 
%     t = (0:length(stim)-1)/fs_stim;
%     subplot(2,2,i)
%     plot(t*1e3,stim(:,i))
%     title(sprintf('Channel %d',i))
% 
% 
% end
% 
% 
% xlabel('Time (ms)')
% ylabel('Amplitude (V)')
% 
% subtitle('Stimulation Channels')

%% Sing looks like the wave to be delivered, with amplitude in uA

%Another attempt to dissect things
%
% timeOfStim = find(Sing1Mask>0);
%
% differenceTimeOfStim = diff(timeOfStim);
% distBetween = find(differenceTimeOfStim>1);
% sampleInd = timeOfStim(distBetween);

%% Sing looks like the wave to be delivered, with amplitude in uA
% Try working from this

% build a burst table with the timing of stimuli
bursts = [];

Sing1 = Sing(:,1);
fs_sing = Sing_info.SamplingRateHz;

samplesOfPulse = round(2*fs_stim/1e3);



% trying something like A_BuildStimTables from BetaStim


Sing1Mask = Sing1~=0;
dmode = diff([0 Sing1Mask' 0 ]);


dmode(end-1) = dmode(end);


bursts(2,:) = find(dmode==1);
bursts(3,:) = find(dmode==-1);

stims = squeeze(getEpochSignal(Sing1,(bursts(2,:)-1),(bursts(3,:))+1));

if strcmp(plotIt,'y')
    t = (0:size(stims,1)-1)/fs_sing;
    t = t*1e3;
    figure
    plot(t,stims)
    xlabel('Time (ms');
    ylabel('Current to be delivered (mA)')
    title('Current to be delivered for all trials')
end

singEpoched = squeeze(getEpochSignal(Sing1,(bursts(2,:)-1),(bursts(3,:))+1));
t = (0:size(singEpoched,1)-1)/fs_sing;
t = t*1e3;

if strcmp(plotIt,'y')
    
    figure
    plot(t,singEpoched)
    xlabel('Time (ms)');
    ylabel('Current to be delivered (\muA)')
    title('Current to be delivered for all trials')
end


% delay loks to be 0.2867 ms from below.

%% Plot stims with info from above

stim1 = stim(:,1);
stim1Epoched = squeeze(getEpochSignal(stim1,(bursts(2,:)-1),(bursts(3,:))+1));
t = (0:size(stim1Epoched,1)-1)/fs_stim;
t = t*1e3;
if strcmp(plotIt,'y')
    figure
    plot(t,stim1Epoched)
    xlabel('Time (ms');
    ylabel('Voltage (V)');
    title('Finding the delay between current output and stim delivery')
    
    hold on
    
    plot(t,stims)
end

% get the delay in stim times

delay = round(0.2867*fs_stim/1e3);

% plot the appropriately delayed signal

if strcmp(plotIt,'y')
    figure
    stimTimesBegin = bursts(2,:)-1+delay;
    stimTimesEnd = bursts(3,:)-1+delay;
    stim1Epoched = squeeze(getEpochSignal(stim1,stimTimesBegin,stimTimesEnd));
    t = (0:size(stim1Epoched,1)-1)/fs_stim;
    t = t*1e3;
    figure
    plot(t,stim1Epoched)
    xlabel('Time (ms');
    ylabel('Voltage (V)');
    title('Stim voltage monitoring with delay added in')
end


%% extract data

% try and account for delay for the stim times
stimTimes = bursts(2,:)-1+delay;

% DJC 7-7-2016, changed presamps and post samps to 1 second
presamps = round(preTime * fs_data / 1e3); % pre time in sec
postsamps = round(postTime * fs_data / 1e3 ); % post time in sec, % modified DJC to look at up to 300 ms after


% sampling rate conversion between stim and data
fac = fs_stim/fs_data;

% find times where stims start in terms of data sampling rate
sts = round(stimTimes / fac);

% looks like there's an additional 14 sample delay between the stimulation being set to
% be delivered....and the ECoG recording. which would be 1.15 ms?

delay2 = 14;
sts = round(stimTimes / fac) + delay2;

%% Interpolation from miah's code

% uncomment this if wanting to interpolate, broken right now

% for i = 1:size(data,2)
%     presamps = round(0.025 * fs_data); % pre time in sec
%     postsamps = round(0.125 * fs_data); % post time in sec, % modified DJC to look at up to 300 ms after
%     eco = data(:,i);
%
%     edd = zeros(size(sts));
%     efs = fs_data;
%
%     temp = squeeze(getEpochSignal(eco, sts-presamps, sts+postsamps+1));
%     foo = mean(temp,2);
%     lastsample = round(0.040 * efs);
%     foo(lastsample:end) = foo(lastsample-1);
%
%     last = find(abs(zscore(foo))>1,1,'last');
%     last2 = find(abs(diff(foo))>30e-6,1,'last')+1;
%
%     zc = false;
%
%     if (isempty(last2))
%         if (isempty(last))
%             error ('something seems wrong in the triggered average');
%         else
%             ct = last;
%         end
%     else
%         if (isempty(last))
%             ct = last2;
%         else
%             ct = max(last, last2);
%         end
%     end
%
%     while (~zc && ct <= length(foo))
%         zc = sign(foo(ct-1)) ~= sign(foo(ct));
%         ct = ct + 1;
%     end
%
%     if (ct > max(last, last2) + 0.10 * efs) % marched along more than 10 msec, probably gone to far
%         ct = max(last, last2);
%     end
%
%     % DJC - 8-31-2015 - i believe this is messing with the resizing
%     % in the figures
%     %             subplot(8,8,chan);
%     %             plot(foo);
%     %             vline(ct);
%     %
%
%     % comment this part out for no interpolation
%     for sti = 1:length(sts)
%         win = (sts(sti)-presamps):(sts(sti)+postsamps+1);
%
%         % interpolation approach
%         eco(win(presamps:(ct-1))) = interp1([presamps-1 ct], eco(win([presamps-1 ct])), presamps:(ct-1));
%     end
%
%     data(:,i) = eco;
% end

%% get the data epochs
dataEpoched = squeeze(getEpochSignal(data,sts-presamps,sts+postsamps+1));

% set the time vector to be set by the pre and post samps
t = (-presamps:postsamps)*1e3/fs_data;

% get average to subtract off
% get average from before signal
if strcmp(sid,'20f8a3')
    epochAve = mean(dataEpoched(t<0,:,:),1);
    meanToSub = repmat(epochAve,size(dataEpoched,1),1,1);
    dataEpoched = dataEpoched -  meanToSub;
end

% if strcmp(sid,'2fd831')
%     epochAve = mean(dataEpoched(t<0,:,:),1);
%     meanToSub = repmat(epochAve,size(dataEpoched,1),1,1);
%     dataEpoched = dataEpoched -  meanToSub;
% end

% ui box for input
prompt = {'sscale the y axis to the maximum stim pulse value? "y" or "n" '};
dlg_title = 'Scale';
num_lines = 1;
defaultans = {'n'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
scaling = answer{1};

if strcmp(scaling,'y')
    maxVal = max(dataEpoched(:));
    minVal = min(dataEpoched(:));
end


%% plot individual trials for each condition on a different graph

labels = max(singEpoched);
uniqueLabels = unique(labels);

% intialize counter for plotting
k = 1;

% make vector of stim channels
stimChans = [stim_1 stim_2];

% determine number of subplot
subPlots = numSubplots(numChans);
p = subPlots(1);
q = subPlots(2);

% plot each condition separately e.g. 1000 uA, 2000 uA, and so on

for i=uniqueLabels
    figure;
    dataInterest = dataEpoched(:,:,labels==i);
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

%% plot averages for 3 conditions on the same graph


k = 1;
figure;
for k = 1:length(dataAvgs)
    
    tempData = dataAvgs{k};
    
    for j = 1:numChans
        s = subplot(p,q,j);
        plot(t,squeeze(tempData(:,j)),'linewidth',2);
        hold on;
        xlim([min(t) max(t)]);
        
        
        % change y axis scaling if necessary
        if strcmp(scaling,'y')
            ylim([minVal maxVal]);
        end
        
        
        if ismember(j,stimChans)
            ax = gca;
            ax.Box = 'on';
            ax.XColor = 'red';
            ax.YColor = 'red';
            ax.LineWidth = 2;
            title(num2str(j),'color','red')
            
        else
            title(num2str(j));
            
        end
        
        vline(0);
        
    end
    gcf;
end
xlabel('time in ms');
ylabel('voltage in V');
subtitle(['Averages for all conditions']);
legLabels = {[num2str(uniqueLabels(1))]};

k = 2;
if length(uniqueLabels>1)
    for i = uniqueLabels(2:end)
        legLabels{end+1} = [num2str(uniqueLabels(k))];
        k = k+1;
    end
end

legend(s,legLabels);

%% 6-23-2016 - plot channel of interest
%
% pick channel
% pick range of stims
j = 1:10;

% if strcmp(plotIt,'y')
%     figure
%     for i = 1:128
%         plot(t,1e6*mean(dataEpoched(:,i,j),3),'m','LineWidth',1)
%         xlabel('time (ms)')
%         ylabel('Amplitude (\muV)')
%         title(['Average for subselected stims for channel ', num2str(i)])
%         ylim([-150 150])
%         xlim([-5 100])
%         pause(1)
%         
%     end
% end
%%
figure
i = 21;

plot(t,1e6*(squeeze(dataEpoched(:,i,j))))
xlabel('time (ms)')
ylabel('Amplitude (\muV)')
title(['Average for subselected stims for channel ', num2str(i)])



%%
%save(fullfile(OUTPUT_DIR, ['stim_constantV',num2str(stim_chans(1)),'_',num2str(stim_chans(2))]), 'data_info','dataEpoched','dataEpochedHigh','dataEpochedLow','dataEpochedMid','fs_data','fs_sing','fs_stim','Sing','Sing_info','stim','stim_chans','stim_info','t');
% save(fullfile(OUTPUT_DIR, ['stim_',num2str(stim_chans(1)),'_',num2str(stim_chans(2))]), 'data_info','dataEpoched','dataEpochedHigh','dataEpochedLow','dataEpochedMid','fs_data','fs_sing','fs_stim','Sing','Sing_info','stim','stim_chans','stim_info','t');
%save(fullfile(OUTPUT_DIR, ['stim_widePulse',num2str(stim_chans(1)),'_',num2str(stim_chans(2))]), 'data_info','dataEpoched','dataEpochedHigh','dataEpochedLow','dataEpochedMid','fs_data','fs_sing','fs_stim','Sing','Sing_info','stim','stim_chans','stim_info','t');
save(fullfile(OUTPUT_DIR, ['stim_widePulse_',num2str(stim_chans(1)),'_',num2str(stim_chans(2))]), 'data_info','dataEpoched','fs_data','fs_sing','fs_stim','Sing','Sing_info','stim','stim_chans','stim_info','t');
