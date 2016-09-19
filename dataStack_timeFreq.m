%% dataStack function for SVD analysis

function [dataStackedGood] = dataStack_timeFreq(dataEpoched,t,post_begin,post_end,channelsOfInt,stim_1,stim_2,bads,fs_data,fw,t_post)
% DATASTACK stacks data for SVD analysis
% This function takes an input structure in frequency x time x channels x trials format
% and stacks all time points, following a beginning time point and before
% an ending one, for a list of good channels, while excluding bad channels
% (always excluding stimulation channels). The next call could for svd or
% something similar

numChans = size(dataEpoched,3);

% use defaults if arguments not supplied in function call
if(~exist('bads','var'))
    bads = [];
end

if(~exist('stim_1','var'))
    stim_1  = [];
end

if(~exist('stim_2','var'))
    stim_2 = [];
end

if(~exist('channelsOfInt','var'))
    channelsOfInt = [1:numChans];
end




% reshift the data so we can do a vector operation to stack all of it like
% we did for that one example

%data_permuted  = permute(dataEpoched,[1,3,2]);
data_permuted  = dataEpoched;

% stack the data

data_stacked_1 = reshape(data_permuted,[size(data_permuted,1)*size(data_permuted,2),size(data_permuted,3),size(data_permuted,4)]);
data_rePerm = permute(data_stacked_1,[1,3,2]);
data_stacked_2 = reshape(data_rePerm,[size(data_rePerm,1)*size(data_rePerm,2),size(data_rePerm,3)]);


% make a vector of all of the channels we have
goods = zeros(numChans,1);

% pick the good channels
goods(channelsOfInt) = 1;

% pick the ones to ignore
badTotal = [stim_1,stim_2,bads];
goods(badTotal) = 0;

% 7-13-2016 - input channels of interest

% make a logical matrix
goods = logical(goods);

% select the good channels
dataStackedGood = data_stacked_2(:,goods);
%% plot example to make sure reshaping worked fine 
% example channel
idx = 60;

% number of trials that were stacked
numTrials = 10;

examp = data_stacked_2(:,idx);
examp_reshape = reshape(examp,[length(fw) length(t_post) numTrials]);

% set colormap using cbrewer
CT = cbrewer('div','RdBu',11);

% flip it so red is increase, blue is down
CT = flipud(CT);

figure
imagesc(t_post, fw, squeeze(examp_reshape(:,:,1))); axis xy
xlabel('Time in ms')
ylabel('Frequency (Hz)')
colorbar
title('Wavelet analysis for Post Stim - Normalized to Pre')
colormap(CT);
set_colormap_threshold(gcf, [-1 1], [-10 10], [.5 .5 .5])

end