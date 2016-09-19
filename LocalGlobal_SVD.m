%% Project local modes onto global modes
close all, clear all, clc

%% load a data file, must update below with the file you want to upload
% Or set allStimPairs to true;
fileName = 'stim_12_52';

allStimPairs = true;
includeLowMed = true; % make true if you want to include the low and medium pulses in addition to the high pulses

% to add paths to all the subfolders
addpath(genpath(pwd))

%SPECIFIC ONLY TO DJC DESKTOP RIGHT NOW
% filePath = 'C:\Users\djcald\Google Drive\GRIDLabDavidShared\CSNE YSP 2016\1sBefore1safter\';
% filePath = 'C:\Users\djcald\GoogleDrive\GRIDLabDavidShared\20f8a3\StimulationSpacingChunked\';

%SPECIFIC ONLY TO JAC DESKTOP RIGHT NOW
filePath = 'C:\Users\jcronin\Data\Subjects\3f2113\data\d6\Matlab\StimulationSpacing\1sBefore1safter\';

%SPECIFIC ONLY TO JAC Laptop RIGHT NOW
% filePath = '/Users/jcronin/Desktop/Data/3f2113/1sBefore1safter/';

if ~allStimPairs
    load([filePath, fileName, '.mat'])
end

%% Start building the data matrices to SVD
% post stim time window
post_begin = 10; % try 10, rather than 5 DJC
post_end = (450+post_begin);

% Build the local and global data matrices

filter_it = 'y';

if ~allStimPairs
    % Disregard bad channels
    prompt = {'what is the list of channels to IGNORE? e.g. 1:8,12 '};
    dlg_title = 'BadChannels';
    num_lines = 1;
    defaultans = {num2str(stim_chans)};
    answerChans = inputdlg(prompt,dlg_title,num_lines,defaultans);
    badChans = str2num(answerChans{1});
    
    goods = ones(size(dataEpochedHigh,2),1);
    goods(badChans) = 0;
    goods(65:end) = 0; % Only look at the grid for now
    goods = logical(goods);
    chansToStack = goods;
    
    if includeLowMed
        X = dataStack(dataEpoched,t,post_begin,post_end,chansToStack,[],[],[],fs_data,filter_it)';
        numPulses = size(dataEpoched, 3);
        stimType = [ones(size(dataEpochedLow,3), 1); 2*ones(size(dataEpochedMid,3), 1); 3*ones(size(dataEpochedMid,3), 1)];
    else
        X = dataStack(dataEpochedHigh,t,post_begin,post_end,chansToStack,[],[],[],fs_data,filter_it)';
        numPulses = size(dataEpochedHigh, 3);
    end
    X_k = reshape(X, [size(X,1), size(X,2)/numPulses, numPulses]);
    X_forSVD = X;
    X_k_forSVD = X_k;
    fileToLoad = 1; % Just need to know for later purposes that you only loaded in one file
else
    % Load in all datafiles one at a time and concatenate them to the X and
    % X_k matrices
    fileNamePrefix = 'stim_';
    stimPairs = {'4_60', '12_52', '20_44' '25_32' '26_31' '27_30' '28_29' '28_36'}; % There's probably a better way to do this, like reading in from the given path...
    fileToLoad = cellfun(@(s) strcat(fileNamePrefix, s), stimPairs, 'UniformOutput', false);
    X = [];
    allStimChans = [];
    stimType = [];
    numPulsesTOTAL = 0;
    for i=1:length(fileToLoad)
        load([filePath, fileToLoad{i}])
        chansToStack = 1:64;
        allStimChans = horzcat(allStimChans, stim_chans);
        if includeLowMed
            X = horzcat(X, dataStack(dataEpoched,t,post_begin,post_end,chansToStack,[],[],[],fs_data,filter_it)');
            numPulsesTOTAL = numPulsesTOTAL + size(dataEpoched, 3);
            stimType = vertcat(stimType, [ones(size(dataEpochedLow,3), 1); 2*ones(size(dataEpochedMid,3), 1); 3*ones(size(dataEpochedMid,3), 1)]);
        else
            X = horzcat(X, dataStack(dataEpochedHigh,t,post_begin,post_end,chansToStack,[],[],[],fs_data,filter_it)');
            numPulsesTOTAL = numPulsesTOTAL + size(dataEpochedHigh, 3);
        end
    end
    X_k = reshape(X, [size(X,1), size(X,2)/numPulsesTOTAL, numPulsesTOTAL]);
    
    % Need to decide how do deal with the fact that we don't want to include
    % the stim channels when they were stimming, but also probably shouldn't
    % discard all 15 stim channels
    
    % First try discarding all of the stim channels:
    goods = ones(size(X,1),1);
    goods(unique(allStimChans)) = 0;
    goods = logical(goods);
    
    X_forSVD = X(goods,:);
    X_k_forSVD = X_k(goods,:,:);
    
    % Or try keeping all channels:
    %     X_forSVD = X;
    %     X_k_forSVD = X_k;
    %     goods = ones(size(X,1),1);
end

%% SVD of the local matrices
% Data is already in m x n form, where m is the number of
% channels, and n is the time points (rows = sensors, columns = samples)
uL=zeros(sum(goods),sum(goods),size(X_k,3));
sL=zeros(size(uL));
vL=zeros(size(X_k,2),sum(goods),size(X_k,3));

for i=1:size(X_k,3)
    [uL(:,:,i), sL(:,:,i), vL(:,:,i)] = svd(X_k_forSVD(:,:,i),'econ');
end

% SVD of the global matrix
[uG, sG, vG] = svd(X_forSVD, 'econ');

modes=1:5;
SVDplot(uG, sG, vG, false, [], modes)


%% Plot projections
modes=[3:5];
a=zeros(1,length(modes));

figure
colorIncrement=0.1;


if includeLowMed
    rcolor=1.0; % this is to control the color of the line
    gcolor=1.0;
    bcolor=1.0;
    for ii=1:size(X_k_forSVD,3)
        ind=1;
        for i=modes
            a(ind) = uG(:,i).'*uL(:,i,ii);
            a_mat(ii,1:3) = a;
            ind=ind+1;
        end
        if stimType(ii)==1 % low stims, red
            if ii~=1 && stimType(ii-1)~=1
                rcolor=1.0; % This resets the color everytime you start a new set of stimType
            end
            scatter3(a(1),a(2),a(3), [], [1.0 0.0 rcolor],'filled')
            rcolor=rcolor-colorIncrement;
            
        elseif stimType(ii)==2 % medium stims, green
            if ii~=1 && stimType(ii-1)~=2
                bcolor=1.0; % This resets the color everytime you start a new set of stimType
            end
            scatter3(a(1),a(2),a(3), [], [bcolor 1.0 0.0],'filled')
            bcolor=bcolor-colorIncrement;
            
        elseif stimType(ii)==3 % high stims, blue
            if ii~=1 && stimType(ii-1)~=3
                gcolor=1.0; % This resets the color everytime you start a new set of stimType
            end
            scatter3(a(1),a(2),a(3), [], [0.0 gcolor 1.0],'filled')
            gcolor=gcolor-colorIncrement;
            
        end
        hold on
    end
    n=1:ceil(sum(stimType==1)/length(fileToLoad));
    labels=arrayfun(@(b) ['Low Pulse ', num2str(b)], n, 'UniformOutput', false);
    n=1:ceil(sum(stimType==2)/length(fileToLoad));
    labels=[labels, arrayfun(@(b) ['Mid Pulse ', num2str(b)], n, 'UniformOutput', false)];
    n=1:ceil(sum(stimType==3)/length(fileToLoad));
    labels=[labels, arrayfun(@(b) ['High Pulse ', num2str(b)], n, 'UniformOutput', false)];
    
else
    labels = cell(1, size(dataEpochedHigh,3));
    for pairs=1:length(fileToLoad)
        gcolor=1.0; % this is to control the color of the line
        for ii=1:size(dataEpochedHigh,3)
            ind=1;
            for i=modes
                a(ind) = uG(:,i).'*uL(:,i,ii+(pairs-1)*size(dataEpochedHigh,3));
                a_mat(pairs,ii,1:3) = a;
                ind=ind+1;
            end
            scatter3(a(1),a(2),a(3), [], [0.0 gcolor 1.0],'filled')
            hold on
            gcolor=gcolor-colorIncrement;
            labels{ii}=['Pulse ', num2str(ii)];
        end
    end
end

xlabel(['Mode ' num2str(modes(1))])
ylabel(['Mode ' num2str(modes(2))])
zlabel(['Mode ' num2str(modes(3))])
legend(labels)

if allStimPairs
    title('All stim pairs: 8 pairs*10 stim pulses, local mode onto global')
else
    title([strrep(fileName, '_', '\_'), ': each point is one stim pulse local mode onto global mode'])
end

