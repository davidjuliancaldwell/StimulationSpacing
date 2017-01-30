function [d_mat,z_mat] = procrustes_metric(t_fullData,t_processedSignal,fullData,processedSignal,pre,post,bads,plotIt)
%USAGE: function [d_mat,z_mat] = procrustes_metric(t_fullData,tNew,fullData,processedSignal,bads,pre,post,plotIt)
%This function will perform a procrustes distance analysis between a signal
%in fullData, and processedSignal. The outputs are d_mat, which is a
% channel by trial matrix showing the goodness of fit metric distance
% between the fullData and the transformed processed signal, and z_mat, which are the
% transformed processedSingal values. 
% NOTE: the two signals can be input with different time vectors, however,
% they both then need to have individual time vectors inserted as shown
% below, as well as both contain the pre and post time values
%
% t_fullData = trial time vector for the full data matrix
% t_processedSignal = trial time vector for the processed signal
% fullData = samples x channels x trials - the full data matrix
% proccesedSignal = the processed signal which will be transformed
% pre = the time point at which to begin extracting the signal
% post = the time point at which to stop extracting the signal
% bads = stimulation channels, or any channels to ignore
% plotIt = plot the procrustes distance results 
% REQUIRES FastICA algorithm in path

% set default plotIt to false
if(~exist('plotIt','var'))
    plotIt = false;
end

% no bad channels unless input 
if(~exist('bads','var'))
    bads = [];
end

% get original data for time vector of interest
data_orig =  fullData((t_fullData>=pre & t_fullData<=post),:,:);

% extract the desired portion of the processed signal
artifact_extract = processedSignal((t_processedSignal>=pre & t_processedSignal<=post),:,:);

% number of trials
numTrials = size(fullData,3);

% get all of the bad channels 
badTotal = bads;

% total channels in data
numChans = size(fullData,2);

% make logical good channels matrix to index
goods = zeros(numChans,1);

channelsOfInt = 1:numChans;

goods(channelsOfInt) = 1;
% set the goods matrix to be zero where bad channels are
goods(badTotal) = 0;
% make it logical
goods = logical(goods);
% this will be used for indexing later 


%%
% initialize matrices 
d_mat = zeros(size(artifact_extract,2),size(numTrials,1));
z_mat = zeros(size(data_orig,1),size(artifact_extract,2),size(numTrials,1));

% channel by channel, trial by trial
for i = 1:numTrials
    
    data_temp = squeeze(data_orig(:,:,i));
    artifact_temp = squeeze(artifact_extract(:,:,i));
    
    for j = 1:size(processedSignal,2)
        
        data_temp_select = data_temp(:,j);
        artifact_temp_select = artifact_temp(:,j);
        
        [d,z,transform] = procrustes(data_temp_select,artifact_temp_select);
        
        d_mat(j,i) = d;
        z_mat(:,j,i) = z;
        
    end
end

d_mat(bads,:) = NaN(length(badTotal),numTrials);
z_mat(:,bads,:) = NaN(length(z), length(badTotal),numTrials);

%% plot procrustes values

if plotIt
    
    % trial plot
    figure
    hold on
    numChans = size(d_mat,1);
    
    for i = 1:numTrials
        
        scatter(repmat(i,size(d_mat,1),1),d_mat(:,i))
        
    end
    xlabel(['Trial Number'])
    ylabel(['Procrustes goodness of fit'])
    title(['Procrustes Distance Metric for Artifact vs. Original Signal'])
    sprintf('Average procrustes GOF for each trial')
    
    average_proc_trial = nanmean(d_mat,1);
    for j = 1:numTrials 
    fprintf('Average procrustes GOF for trial # %d = %0.4f\n',j,average_proc_trial(j))
    end
    %
    figure
    hold on
    for i = 1:numChans
        
        scatter(repmat(i,size(d_mat,2),1),d_mat(i,:))
        
    end
    xlabel(['Channel Number'])
    ylabel(['Procrustes goodness of fit'])
    title(['Procrustes Distance Metric for Artifact vs. Original Signal'])
    
    average_proc_channel = nanmean(d_mat,2);
    for j = 1:numChans
    fprintf('Average procrustes GOF for channel # %d = %0.4f\n',j,average_proc_channel(j))
    end
end

end