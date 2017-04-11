function [epochedSignal, t] = stimulationExtractionBetween(data, fs_data, stimTrig, fs_stim, L)
%% 3-31-2017 - Jeneva Cronin - script to extract data between stims
% data that you wish to extract from
% stimTrig: the stimulation signal that triggered the stim (usually this will
% be the Sing.data signal
% L: desired length of epochs
% Returns:
% epochedSignal - the signal between stim trains, length of L
% t - the timing matrix; each row corresponds to a specific epoch and is
% the time in ms after the end of the burst


plotIt = 'y';

% build a burst table with the timing of stimuli
bursts = [];

% Sing1 = Sing.data(:,1);
% fs_sing = Sing_info.SamplingRateHz;

% samplesOfPulse = round(2*stim_fs/1e3);


% trying something like A_BuildStimTables from BetaStim
stimMask = stim~=0;
dmode = diff([0 stimMask' 0 ]);
dmode(end-1) = dmode(end);

bursts(2,:) = find(dmode==1); % When the stim pulse starts
bursts(3,:) = find(dmode==-1); % When the stim pulse ends

stims = squeeze(getEpochSignal(stimTrig,(bursts(2,:)-1),(bursts(3,:))+1));

if strcmp(plotIt,'y')
    t = (0:size(stims,1)-1)/fs_sing;
    t = t*1e3;
    figure
    plot(t,stims)
    xlabel('Time (ms');
    ylabel('Current to be delivered (mA)')
    title('Current to be delivered for all trials')
end

%% Determine the middle time between stims
% This subtracts the previous ending burst sample number from the starting
% burst sample number, to get the samples between the bursts. Then divide
% by 2 and add back to the ending burst sample to get the middle sample. 
sampsBetween = (bursts(2,2:end) - bursts(3,1:end-1));
bursts(4,2:end) = floor(sampsBetween/2) + bursts(3,1:end-1); 
% Set the first 'middle sample' to be mode(sampsBetween/2) before the first
% burst 
bursts(4,1) = bursts(2,1) - mode(floor(sampsBetween/2));

%% Determine lengths of epochs
% Decide how much data to pull out before and after the middle sample
% value. This will be equal to L/2 (so that the whole epoch length is L),
% unless that would cause it to include stimulations, in which case, L must
% be shortened
if L >= min(floor(sampsBetween))
    % L is too long, so need to shorten
    L = min(floor(sampsBetween));
    warning('Shortened the length of the epochs (parameter L)')
end 

%% Set stim delivery delay based on previous work:
% get the delay in stim times (delay is about 0.2867 ms)
delay = round(0.2867*fs_stim/1e3) + 14;

%% Extract data
% try and account for delay for the stim times
betweenTimes = bursts(4,:)-1+delay;

% sampling rate conversion between stim and data
fac = fs_stim/fs_data;

% find times between stimulation samples in terms of data sampling rate
% (right now we're in terms of stim sampling rate)
sts = round(betweenTimes / fac);

epochedSignal = squeeze(getEpochSignal(data,sts-L/2,sts+L/2));

%% Set the time matrix based on post stim timing
t = (-L/2:postsamps)*1e3/fs_data;
floor(sampsBetween/2)

% get average to subtract off

if strcmp(sid,'20f8a3')
    epochAve = mean(dataEpoched,1);
    meanToSub = repmat(epochAve,size(dataEpoched,1),1,1);
    dataEpoched = dataEpoched -  meanToSub;
end



%% 6-23-2016
% plot aggregate data for each stim type - THIS ONLY WORKS IF THERE ARE 30
% TOTAL STIM EPOCHS - stim 9 for instance only has 28 total

% chunk out data

% to separate out low and high

%k=1:8;
%j = 1:8;

k = 1:10;
j = 1:10;

% figure
for i = 1:64
    %     hold on
    %     subplot(8,8,i)
    %     plot(t,1e6*mean(dataEpoched(:,i,j),3),'m','LineWidth',1)
    %
    %     ylim([-150 150])
    %     xlim([-100 300])
    %     title(sprintf('Channel %d',i))
    %     %     pause(1)
    
    dataEpochedLow(:,i,k) = dataEpoched(:,i,j);
    
end
% subtitle('Average traces for all stimulations - means not subtracted - stims 1:10')
% xlabel('Time (ms)')
% ylabel('Voltage (uV)')
%
% figure
for i = 65:128
    %     hold on
    %     subplot(8,8,i-64)
    %     plot(t,1e6*mean(dataEpoched(:,i,j),3),'m','LineWidth',1)
    %
    %     ylim([-150 150])
    %     xlim([-100 300])
    %     title(sprintf('Channel %d',i))
    %     %     pause(1)
    dataEpochedLow(:,i,k) = dataEpoched(:,i,j);
    
    
end
% subtitle('Average traces for all stimulations - means not subtracted - stims 1:10')
% xlabel('Time (ms)')
% ylabel('Voltage (uV)')

k= 1:10;
%j = 9:18;
j = 11:20;

% figure
for i = 1:64
    %     hold on
    %     subplot(8,8,i)
    %     plot(t,1e6*mean(dataEpoched(:,i,j),3),'m','LineWidth',1)
    %
    %     ylim([-150 150])
    %     xlim([-100 300])
    %     title(sprintf('Channel %d',i))
    %     %     pause(1)
    dataEpochedMid(:,i,k) = dataEpoched(:,i,j);
    
    
end
% subtitle('Average traces for all stimulations - means not subtracted stims - stims 11:20')
% xlabel('Time (ms)')
% ylabel('Amplitude (\muV)')
%
% figure
for i = 65:128
    %     hold on
    %     subplot(8,8,i-64)
    %     plot(t,1e6*mean(dataEpoched(:,i,j),3),'m','LineWidth',1)
    %
    %     ylim([-150 150])
    %     xlim([-100 300])
    %     title(sprintf('Channel %d',i))
    %     %     pause(1)
    dataEpochedMid(:,i,k) = dataEpoched(:,i,j);
    
    
end
% subtitle('Average traces for all stimulations - means not subtracted - stims 11:20 ')
% xlabel('Time (ms)')
% ylabel('Amplitude (\muV)')

%k=1:10;
j = 21:37;
%j = 21:39;

% figure
for i = 1:64
    %     hold on
    %     subplot(8,8,i)
    %     plot(t,1e6*mean(dataEpoched(:,i,j),3),'m','LineWidth',1)
    %
    %     ylim([-150 150])
    %     xlim([-100 300])
    %     title(sprintf('Channel %d',i))
    %     %     pause(1)
    dataEpochedHigh(:,i,k) = dataEpoched(:,i,j);
    
    
end
% subtitle('Average traces for all stimulations - means not subtracted - stims 21:30 ' )
% xlabel('Time (ms)')
% ylabel('Amplitude (\muV)')
%
% figure
for i = 65:128
    %     hold on
    %     subplot(8,8,i-64)
    %     plot(t,1e6*mean(dataEpoched(:,i,j),3),'m','LineWidth',1)
    %
    %     ylim([-150 150])
    %     xlim([-100 300])
    %     title(sprintf('Channel %d',i))
    %     %     pause(1)
    dataEpochedHigh(:,i,k) = dataEpoched(:,i,j);
    
    
end
% subtitle('Average traces for all stimulations - means not subtracted - stims 21:30')
% xlabel('Time (ms)')
% ylabel('Amplitude (\muV)')


%% 6-23-2016 - plot channel of interest
%
% pick channel
i = 21;
% pick range of stims
%j = 1:10;
%j = 11:20;
j = 21:30;

figure
for i = 1:128
    plot(t,1e6*mean(dataEpoched(:,i,j),3),'m','LineWidth',1)
    xlabel('time (ms)')
    ylabel('Amplitude (\muV)')
    title(['Average for subselected stims for channel ', num2str(i)])
    ylim([-150 150])
    xlim([-5 100])
    pause(1)
    
end
% figure
% plot(t,1e6*(squeeze(dataEpoched(:,i,j))))
% xlabel('time (ms)')
% ylabel('Amplitude (\muV)')
% title(['Average for subselected stims for channel ', num2str(i)])


end