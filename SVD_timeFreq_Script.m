%% script to do SVD on time-frequency plots

close all, clear all, clc

%% load a data file, must update below with the file you want to upload
% Or set allStimPairs to true;
fileName = 'stim_12_52';

allStimPairs = false;
includeLowMed = false; % make true if you want to include the low and medium pulses in addition to the high pulses

% to add paths to all the subfolders
addpath(genpath(pwd))

%SPECIFIC ONLY TO DJC DESKTOP RIGHT NOW
filePath = 'C:\Users\djcald\Google Drive\GRIDLabDavidShared\CSNE YSP 2016\1sBefore1safter\';
% filePath = 'C:\Users\djcald\GoogleDrive\GRIDLabDavidShared\20f8a3\StimulationSpacingChunked\';

%SPECIFIC ONLY TO JAC DESKTOP RIGHT NOW
%filePath = 'C:\Users\jcronin\Data\Subjects\3f2113\data\d6\Matlab\StimulationSpacing\1sBefore1safter\';

%SPECIFIC ONLY TO JAC Laptop RIGHT NOW
% filePath = '/Users/jcronin/Desktop/Data/3f2113/1sBefore1safter/';


%% load it in!

if ~allStimPairs
    load([filePath, fileName, '.mat'])
    
    prompt = {'what is the list of channels to IGNORE? e.g. 1:8,12 ','use sig CCEPs?'};
    dlg_title = 'BadChannels';
    num_lines = 1;
    defaultans = {num2str(stim_chans),'n'};
    answerChans = inputdlg(prompt,dlg_title,num_lines,defaultans);
    badChans = str2num(answerChans{1});
    sigCCEPusage = answerChans{2};
    
    goods = ones(size(dataEpochedHigh,2),1);
    goods(badChans) = 0;
    goods(65:end) = 0; % Only look at the grid for now
    goods = logical(goods);
    chansToStack = goods;
    
    %% initialize parameters
    
    % pre time window
    % try making smaller for wavelet 9-18-2016
    pre_begin = -250;
    pre_end = 0;
    % post time window
    post_begin = 8;
    post_end = (250+post_begin);
    
    % extract pre
    t_pre = t(t<pre_end & t>pre_begin);
    
    % extract post
    t_post = t(t>post_begin & t<post_end);
    
    filter_it = 'y';
    plotIt = false;
    
    %%
    % z-score threshold
    % ui box for input - pick zscore threshold
    
    if strcmp(sigCCEPusage,'y')
        % pre time window
        pre_beginSig = -450;
        pre_endSig = 0;
        % post time window
        post_beginSig = 5;
        post_endSig = (450+post_beginSig);
        
        sigCCEPs = sigCCEP_extract(dataEpochedHigh,t,pre_beginSig,pre_endSig,post_beginSig,post_endSig,fs_data,filter_it,[]);
    end
    
    
    %% wavelets
    signal = dataEpochedHigh;
    [ncPost_m,t_post,fw] = make_waveletMatrix(signal,t,pre_end,pre_begin,post_begin,post_end,fs_data,filter_it,[]);
    
    
    %% subselect channels to stack
    
    if strcmp(sigCCEPusage,'y')
        chansToStack = zeros(size(ncPost_m,3),1);
        chansToStack(sigCCEPs) = 1;
        chansToStack(65:end) = 0;
        chansToStack = logical(chansToStack);
    else
        goods = ones(size(ncPost_m,3),1);
        goods(badChans) = 0;
        goods(65:end) = 0; % Only look at the grid for now
        goods = logical(goods);
        chansToStack = goods;
    end
    
    
    %% stack it
    dataStacked_timeFreq = dataStack_timeFreq(ncPost_m,t,post_begin,post_end,chansToStack,[],[],badChans,fs_data,fw,t_post);
    
    %% svd
    [u,s,v] = svd(dataStacked_timeFreq','econ');
    
    % optimal singular value
    
    m = size(dataStacked_timeFreq',1);
    n = size(dataStacked_timeFreq',2);
    
    % this y is just the singular values sorted
    
    diag(s);
    
    % this is the singular value cutoff
    coef_S_truncation = optimal_SVHT_coef(m/n,0)*median(diag(s));
    
    %%
    figure
    plot(diag(s),'ko','Linewidth',[2])
    % to get percentage in mode
    subplot(2,1,1) % plot normal
    plot(diag(s)/sum(diag(s)),'ko','Linewidth',[2])
    title('singular values, fractions')
    set(gca,'fontsize',14)
    
    subplot(2,1,2) % plot semilog
    semilogy(diag(s)/sum(diag(s)),'ko','Linewidth',[2])
    title('singular values, fractions, semilog plot')
    set(gca,'fontsize',14)
    
    if exist('coef_S_truncation', 'var')
        figure
        
        plot(diag(s),'ko','Linewidth',[2])
        hold on
        
        plot([0 length(s)], [coef_S_truncation coef_S_truncation], 'r')
        title('singular values and cutoff')
        
    end
    %%
    % reshape and plot modes
    % set colormap using cbrewer
    CT = cbrewer('div','RdBu',11);
    
    % flip it so red is increase, blue is down
    CT = flipud(CT);
    figure
    modes = [1:3];
    numSubs = length(modes);
    for i=1:length(modes)
        subplot(numSubs,1,i)
        imagesc(t_post,fw,(reshape(v(1:length(v)/10,modes(i)),[length(fw) length(t_post)]))); axis xy
        title(['Temporal portion of mode #: ', num2str(modes(i))])
        xlabel('Time in ms')
        ylabel('Frequency (Hz)')
        colorbar
        colormap(CT);
    end
    
    %% plot all ten
    
    figure
    
    t = 1e3*[0:length(v)-1]/fs_data;
    
    for i=1:length(modes)
        subplot(numSubs,1,i)
        v_mode = v(:,modes(i));
        v_mode_reshape = reshape(v_mode,[length(v_mode)/size(signal,3),size(signal,3)]);
        v_sequential = reshape(v_mode_reshape,[length(fw),length(t_post)*size(signal,3),1]);
        imagesc(t,fw,(v_sequential)); axis xy
        title(['Temporal portion of mode #: ', num2str(modes(i))])
        xlabel('Time in ms')
        ylabel('Frequency (Hz)')
        colorbar
        colormap(CT);
    end
    
    
    X = dataStack(dataEpochedHigh,t,post_begin,post_end,chansToStack,[],[],[],fs_data,filter_it)';
    
    
    X_k = reshape(X, [size(X,1), size(X,2)/size(dataEpochedHigh, 3), size(dataEpochedHigh, 3)]);
    X_forSVD = X;
    X_k_forSVD = X_k;
    fileToLoad = 1; % Just need to know for later purposes that you only loaded in one file
    
    
    % if ~allStimPairs
    %     % Disregard bad channels
    %     prompt = {'what is the list of channels to IGNORE? e.g. 1:8,12 '};
    %     dlg_title = 'BadChannels';
    %     num_lines = 1;
    %     defaultans = {num2str(stim_chans)};
    %     answerChans = inputdlg(prompt,dlg_title,num_lines,defaultans);
    %     badChans = str2num(answerChans{1});
    %
    %     goods = ones(size(dataEpochedHigh,2),1);
    %     goods(badChans) = 0;
    %     goods(65:end) = 0; % Only look at the grid for now
    %     goods = logical(goods);
    %     chansToStack = goods;
    %
    %     if includeLowMed
    %         X = dataStack(dataEpoched,t,post_begin,post_end,chansToStack,[],[],[],fs_data,filter_it)';
    %         numPulses = size(dataEpoched, 3);
    %         stimType = [ones(size(dataEpochedLow,3), 1); 2*ones(size(dataEpochedMid,3), 1); 3*ones(size(dataEpochedMid,3), 1)];
    %     else
    %         X = dataStack(dataEpochedHigh,t,post_begin,post_end,chansToStack,[],[],[],fs_data,filter_it)';
    %         numPulses = size(dataEpochedHigh, 3);
    %     end
    %     X_k = reshape(X, [size(X,1), size(X,2)/numPulses, numPulses]);
    %     X_forSVD = X;
    %     X_k_forSVD = X_k;
    %     fileToLoad = 1; % Just need to know for later purposes that you only loaded in one file
else
    % Load in all datafiles one at a time and concatenate them to the X and
    % X_k matrices
    fileNamePrefix = 'stim_';
    stimPairs = {'4_60', '12_52', '20_44' '25_32' '26_31' '27_30' '28_29' '28_36'}; % There's probably a better way to do this, like reading in from the given path...
    fileToLoad = cellfun(@(s) strcat(fileNamePrefix, s), stimPairs, 'UniformOutput', false);
    
    %% initialize parameters
    
    % pre time window
    % try making smaller for wavelet 9-18-2016
    pre_begin = -250;
    pre_end = 0;
    % post time window
    post_begin = 8;
    post_end = (250+post_begin);
    
    % extract pre
    t_pre = t(t<pre_end & t>pre_begin);
    
    % extract post
    t_post = t(t>post_begin & t<post_end);
    
    filter_it = 'y';
    plotIt = false;
    
    
    X = [];
    allStimChans = [];
    stimType = [];
    numPulsesTOTAL = 0;
    for i=1:length(fileToLoad)
        load([filePath, fileToLoad{i}])
        chansToStack = 1:64;
        allStimChans = horzcat(allStimChans, stim_chans);
        X = horzcat(X, dataStack(dataEpochedHigh,t,post_begin,post_end,chansToStack,[],[],[],fs_data,filter_it)');
        numPulsesTOTAL = numPulsesTOTAL + size(dataEpochedHigh, 3);
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
