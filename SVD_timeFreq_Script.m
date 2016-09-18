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
pre_begin = -450;
pre_end = 0;
% post time window
post_begin = 5;
post_end = (450+post_begin);

% extract pre
t_pre = t(t<pre_end & t>pre_begin);

% extract post
t_post = t(t>post_begin & t<post_end);

filter_it = 'y';
plotIt = false;

signal = dataEpochedHigh(:,goods,:);

% pre initialize 
chans = numChans(goods);

ncPost_m = zeros(200,length(t_post),length(chans),size(dataEpochedHigh,3));

for i = numChans(goods)
    
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
        ncPost_m(:,:,j,i) = ncPost;
        
    end
end
