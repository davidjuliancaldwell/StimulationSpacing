%% procrustes analysis as metric of how similar stimulus is

% d = procrustes(X,Y) determines a linear transformation (translation, reflection, orthogonal rotation, and scaling) of the points in matrix Y to best conform them to the points in matrix X. The goodness-of-fit criterion is the sum of squared errors. procrustes returns the minimized value of this dissimilarity measure in d. d is standardized by a measure of the scale of X, given by:
%
% sum(sum((X-repmat(mean(X,1),size(X,1),1)).^2,1))
% That is, the sum of squared elements of a centered version of X. However, if X comprises repetitions of the same point, the sum of squared errors is not standardized.
%
% X and Y must have the same number of points (rows), and procrustes matches Y(i) to X(i). Points in Y can have smaller dimension (number of columns) than those in X. In this case, procrustes adds columns of zeros to Y as necessary.
%
% [d,Z] = procrustes(X,Y) also returns the transformed Y values.
%
% [d,Z,transform] = procrustes(X,Y) also returns the transformation that maps Y to Z. transform is a structure array with fields:
%
% c — Translation component
% T — Orthogonal rotation and reflection component
% b — Scale component
% That is:
%
% c = transform.c;
% T = transform.T;
% b = transform.b;
%
% Z = b*Y*T + c;
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% function to call
%[d,Z,transform] = procrustes(X,Y); % also returns the transformation that maps Y to Z. transform is a structure array with fields:

%clear work space
close all;clear all;clc

%%
% load data
% DJC desktop

load('C:\Users\djcald\Google Drive\GRIDLabDavidShared\CSNE YSP 2016\1sBefore1safter\stim_12_52.mat')
%load('C:\Users\djcald\Google Drive\GRIDLabDavidShared\CSNE YSP 2016\1sBefore1safter\stim_28_29.mat')


%%
% extract data of interest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pre time window
pre_begin = -400;
pre_end = 0;
% post time window
post_begin = 5;

% used 100, then tried 20
post_end = (450+post_begin);
% extract pre
t_pre = t(t<pre_end & t>pre_begin);

% extract post
t_post = t(t>post_begin & t<post_end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get stim channels
stimChan1 = stim_chans(1);
stimChan2= stim_chans(2);

%
% make bad channels vector
bads = [];
badTotal = [stimChan1, stimChan2, bads];

% total channels
numChans = 64;
% make logical good channels matrix to index
goods = zeros(numChans,1);

channelsOfInt = 1:numChans;

goods(channelsOfInt) = 1;
% set the goods matrix to be zero where bad channels are
goods(badTotal) = 0;
% make it logical
goods = logical(goods);

% extract the data of interest
dataInt = dataEpochedHigh(:,goods,:);
% NOTE THIS IS DIFFERENT THAN BEFORE, WE WANT TO KEEP STIMULATION IN THERE
dataIntTime = dataInt((t>pre_begin & t<post_end),:,:);

timeVec = t(t>pre_begin & t<post_end);
numTrials = size(dataIntTime,3);

data_permuted  = permute(dataIntTime,[1,3,2]);

% stack the data

data_stacked = reshape(data_permuted,[size(data_permuted,1)*size(data_permuted,2),size(data_permuted,3)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%channelInt = 21;
channelInt = 58;

figure
plot(data_stacked(:,channelInt))
title(['Original Data for Channel ', num2str(channelInt)])

%% Try single trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_icasig = {};
i_mixing_mat = {};
i_sep_mat = {};
i_icasigS = {};
i_mixing_matS = {};
i_sep_matS = {};

scale_factor = 1000;
% visualize trial by trial across grid
%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:size(dataIntTime,3)
    sig_epoch = squeeze(dataIntTime(:,:,i));
    figure
    plot(sig_epoch)
    title(['Trial ' num2str(i)])
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ica across
for i = 1:size(dataIntTime,3)
    sig_epoch = scale_factor.*squeeze(dataIntTime(:,:,i));
    [icasig_temp,mixing_mat_temp,sep_mat_temp] = fastica(sig_epoch');
    
    i_icasigS{i} = icasig_temp;
    i_mixing_matS{i} = mixing_mat_temp;
    i_sep_matS{i} = sep_mat_temp;
    
end

numInt = min(size(icasig_temp,1),10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% visualize ICA
for j = 1:size(dataIntTime,3)
    figure
    for i = 1:numInt
        subplot(numInt,1,i)
        plot(i_icasigS{j}(i,:))
        title(['ICA component # , ', num2str(i)])
    end
    subtitle(['Trial # ', num2str(j)])
    
end

recon_artifact = {};

%%%%%%%%%%%%%%%%%%%%%%%
% reconstruct stim artifact across channels
% try 8 for 1000, 6 for 100 scale factor, 4 for 10 scale factor
%%
num_modes_subtract = 10;
recon_artifact_matrix = [];

%%

% make matrix of reconstruction artifacts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:size(dataIntTime,3)
    
    
    recon_artifact_temp = (i_mixing_matS{i}(:,1:num_modes_subtract)*i_icasigS{i}(1:num_modes_subtract,:))'./scale_factor;
    
    recon_artifact{i} = recon_artifact_temp;
    recon_artifact_matrix(:,:,i) = recon_artifact_temp;
    
    figure
    plot(recon_artifact_temp(:,channelInt))
    hold on
    plot(dataIntTime(:,channelInt,i))
    title(['Channel ', num2str(channelInt), ' Trial ', num2str(i), 'Number of ICA modes subtracted = ', num2str(num_modes_subtract)])
    legend({'recon artifact','original signal'})
    
end


%% this is from ICA prototype
%% single trial scaled - BEST THUS FAR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


i_icasigS = {};
i_mixing_matS = {};
i_sep_matS = {};

scale_factor = 1000;


for i = 1:size(dataIntTime,3)
    sig_epoch = scale_factor.*squeeze(dataIntTime(:,:,i));
    [icasig_temp,mixing_mat_temp,sep_mat_temp] = fastica(sig_epoch');
    
    i_icasigS{i} = icasig_temp;
    i_mixing_matS{i} = mixing_mat_temp;
    i_sep_matS{i} = sep_mat_temp;
    
end


%% visualize the trial by trial ICA components

%
% numInt = min(size(icasig_temp,1),10);
%
% for j = 1:size(dataIntTime,3)
%     figure
%     for i = 1:numInt
%         subplot(numInt,1,i)
%         plot(i_icasigS{j}(i,:))
%         title(['ICA component # , ', num2str(i)])
%     end
%     subtitle(['Trial # ', num2str(j)])
%
% end

%% subtract each one of these components

subtracted_sig_cellS_I = {};

% try 8 for 1000, 6 for 100 scale factor, 4 for 10 scale factor
num_modes_subtract = 10;
subtracted_sig_matrixS_I = [];

for i = 1:size(dataIntTime,3)
    
    
    combined_ica_recon = (i_mixing_matS{i}(:,1:num_modes_subtract)*i_icasigS{i}(1:num_modes_subtract,:))';
    
    % subtracted_sig_ICA_temp = dataIntTime(:,:,i) - combined_ica_recon./scale_factor;
    subtracted_sig_ICA_temp = dataIntTime(:,:,i) - combined_ica_recon./scale_factor;
    
    subtracted_sig_cellS_I{i} = subtracted_sig_ICA_temp;
    subtracted_sig_matrixS_I(:,:,i) = subtracted_sig_ICA_temp;
    
    %     figure
    %     plot(subtracted_sig_ICA_temp(:,channelInt))
    %     hold on
    %     plot(dataIntTime(:,channelInt,i))
    %     title(['Channel ', num2str(channelInt), ' Trial ', num2str(i), 'Number of ICA modes subtracted = ', num2str(num_modes_subtract)])
    %     legend({'subtracted signal','original signal'})
    %
    %     figure
    %     plot(subtracted_sig_ICA_temp(:,channelInt))
    %     title(['Subtracted Signal for ', num2str(num_modes_subtract), ' ICA modes, Channel ', num2str(channelInt), ' Trial ', num2str(i)])
    %
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% procrustes

pre_procrust = -5;
post_procrust = 40;

data_orig =  dataIntTime((timeVec>pre_procrust & timeVec<post_procrust),:,:);
artifact = recon_artifact_matrix((timeVec>pre_procrust & timeVec<post_procrust),:,:);
procrustes_mat = zeros(size(artifact));

%[d,Z,transform] = procrustes(X,Y); % also returns the transformation that maps Y to Z. transform is a structure array with fields:
%%
% channel by channel, trial by trial
for i = 1:size(artifact,3)
    
    data_temp = squeeze(data_orig(:,:,i));
    artifact_temp = squeeze(artifact(:,:,i));
    
    
    
    for j = 1:size(artifact,2)
        
        data_temp_select = data_temp(:,j);
        artifact_temp_select = artifact_temp(:,j);
        
        [d,Z,transform] = procrustes(data_temp_select,artifact_temp_select);
        
        d_mat(j,i) = d;
        Z_mat(:,j,i) = Z;
        
        % transform_cell{j}{i} = transform;
        
        
    end
    
end


%% trial by trial, all channels at once procrustes

for i = 1:size(artifact,3)
    
    data_temp = squeeze(data_orig(:,:,i));
    artifact_temp = squeeze(artifact(:,:,i));
    
    [d,Z,transform] = procrustes(data_temp,artifact_temp);
    
    d_mat_allChans(i) = d;
    Z_mat_allChans(:,:,i) = Z;
    
    transform_cell{i} = transform;
    
    
    
end

