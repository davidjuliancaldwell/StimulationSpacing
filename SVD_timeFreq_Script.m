%% script to do SVD on time-frequency plots

close all, clear all, clc

%% load a data file, must update below with the file you want to upload
% Or set allStimPairs to true;
fileName = 'stim_12_52';

allStimPairs = false;

% to add paths to all the subfolders
addpath(genpath(pwd))

%SPECIFIC ONLY TO DJC DESKTOP RIGHT NOW
filePath = 'C:\Users\djcald\Google Drive\GRIDLabDavidShared\CSNE YSP 2016\1sBefore1safter\';
% filePath = 'C:\Users\djcald\GoogleDrive\GRIDLabDavidShared\20f8a3\StimulationSpacingChunked\';

%SPECIFIC ONLY TO JAC DESKTOP RIGHT NOW
%filePath = 'C:\Users\jcronin\Data\Subjects\3f2113\data\d6\Matlab\StimulationSpacing\1sBefore1safter\';

%SPECIFIC ONLY TO JAC Laptop RIGHT NOW
% filePath = '/Users/jcronin/Desktop/Data/3f2113/1sBefore1safter/';

if ~allStimPairs
    load([filePath, fileName, '.mat'])
end

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
numChans = 1:size(dataEpoched,2);
%     X = dataStack(dataEpochedHigh,t,post_begin,post_end,chansToStack,[],[],[],fs_data,filter_it)';
%     X_k = reshape(X, [size(X,1), size(X,2)/size(dataEpochedHigh, 3), size(dataEpochedHigh, 3)]);
%     X_forSVD = X;
%     X_k_forSVD = X_k;
%     fileToLoad = 1; % Just need to know for later purposes that you only loaded in one file
%
%%

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

signal = dataEpochedHigh;

% pre initialize
chans = numChans(goods);
freqs = [1:3:200];
numFreqs = length(freqs);


ncPost_m = zeros(numFreqs,length(t_post),length(chans),size(dataEpochedHigh,3));

for i = chans
    
    sig = squeeze(signal(:,i,:));
    
    for j = 1:size(sig,2)
        
        if strcmp(filter_it,'y')
            sig_pre = notch(sig((t<pre_end & t>pre_begin),j),[60 120 180 240],fs_data);
            sig_postL = notch(sig((t>post_begin & t<post_end),j),[60 120 180 240],fs_data);
        else
            
            sig_pre = sig((t<pre_end & t>pre_begin),j);
            sig_postL = (sig((t>post_begin & t<post_end),j));
            
        end
        
        if plotIt == true
            figure
            [f_pre,P1_pre] = spectralAnalysis(fs_data,t_pre,sig_pre);
            [f_postL,P1_postL] = spectralAnalysis(fs_data,t_post,sig_postL);
            
            legend({'pre','high'})
        end
        % do some time frequency analysis
        
        if plotIt == true
            figure
        end
        
        [t_post,fw,ncPost] = timeFrequencyAnalWavelet(sig_pre,sig_postL,t_pre,t_post,fs_data,plotIt);
        ncPost_m(:,:,i,j) = ncPost;
        
    end
end

%% once have ncPost_m, stack it, SVD it, etc

goods = ones(size(ncPost_m,3),1);
goods(badChans) = 0;
goods(65:end) = 0; % Only look at the grid for now
goods = logical(goods);
chansToStack = goods;

dataStacked_timeFreq = dataStack_timeFreq(ncPost_m,t,post_begin,post_end,chansToStack,[],[],badChans,fs_data,fw,t_post);

[u,s,v] = svd(dataStacked_timeFreq','econ');

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


