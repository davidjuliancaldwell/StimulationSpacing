%% DJC -12-12-2016
% prototyping with ICA for stimulation extraction prior to data analysis
%%

% clear work space
close all;clear all;clc

%%
% load data
% DJC desktop

load('C:\Users\djcald\Google Drive\GRIDLabDavidShared\CSNE YSP 2016\1sBefore1safter\stim_12_52.mat')

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

numTrials = size(dataIntTime,3);

data_permuted  = permute(dataIntTime,[1,3,2]);

% stack the data

data_stacked = reshape(data_permuted,[size(data_permuted,1)*size(data_permuted,2),size(data_permuted,3)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

channelInt = 58;

figure
plot(data_stacked(:,channelInt))
title(['Original Data for Channel ', num2str(channelInt)])


%% ICA THAT !
% stacked data
%%%%%%%%%%%%%%%%%%%%%
[icasig,mixing_mat,sep_mat] = fastica(data_stacked');


% subtract however many ICA components

% 3 for 25 ms time after, 2 for 100 ms after
num_modes_subtract = 2;

numInt = min(size(icasig,1),10);

figure
for i = 1:numInt
    subplot(numInt,1,i)
    plot(icasig(i,:))
    title(['ICA component # , ', num2str(i)])
end


combined_ica_recon = (mixing_mat(:,1:num_modes_subtract)*icasig(1:num_modes_subtract,:))';

subtracted_sig_ICA = data_stacked - combined_ica_recon;

figure
plot(subtracted_sig_ICA(:,channelInt))
hold on
plot(data_stacked(:,channelInt))
title(['Channel ' num2str(channelInt), 'Number of ICA modes subtracted = ', num2str(num_modes_subtract)])
legend({'subtracted signal','original signal'})

figure
plot(subtracted_sig_ICA(:,channelInt))
title(['Subtracted Signal for ', num2str(num_modes_subtract), ' ICA modes, Channel ', num2str(channelInt)])

%% add on rPCA

% the corrupt data matrix is decomposed into real and imaginary parts
% before applying the inexact_alm_rpca algorithm
scale_fac = 1;

ur=real(scale_fac.*data_stacked.');
ui=imag(scale_fac.*data_stacked.');

% The parameter lambda used in the inexact_alm_rpca algorithm can be tuned
% to best separate the sparse from low-rank matrices
lambda=0.7; 
% See below for examples of changing lambda

% low-rank (R1) and sparse (R2) portions
[R1r, R2r] = inexact_alm_rpca(real(ur.'),lambda);
[R1i, R2i] = inexact_alm_rpca(real(ui.'),lambda);

% the real and imaginary portions are put back together here to SVD the
% low-rank portion:
R1 = R1r + 1i*R1i;
R2 = R2r + 1i*R2i;

[U3,S3,V3]=svd(R1.');



%% try scaled ICA - for all trials stacked

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scale_fac = 1000;
[icasigS,mixing_matS,sep_matS] = fastica(scale_fac.*data_stacked');


numInt = min(size(icasigS,1),10);

%%
num_modes_subtract = 8;

figure
for i = 1:numInt
    subplot(numInt,1,i)
    plot(icasigS(i,:))
    title(['ICA component # , ', num2str(i)])
end


combined_ica_reconS = (mixing_matS(:,1:num_modes_subtract)*icasigS(1:num_modes_subtract,:))';

subtracted_sigS_ICA = data_stacked - combined_ica_reconS./scale_fac;

figure
plot(subtracted_sigS_ICA(:,channelInt))
hold on
plot(data_stacked(:,channelInt))
title(['Channel ' num2str(channelInt), 'Number of ICA modes subtracted = ', num2str(num_modes_subtract)])
legend({'subtracted signal','original signal'})

figure
plot(subtracted_sigS_ICA(:,channelInt))
title(['Subtracted Signal for ', num2str(num_modes_subtract), ' ICA modes, Channel ', num2str(channelInt)])

%% Try single trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_icasig = {};
i_mixing_mat = {};
i_sep_mat = {};


for i = 1:size(dataIntTime,3)
    sig_epoch = squeeze(dataIntTime(:,:,i));
    [icasig_temp,mixing_mat_temp,sep_mat_temp] = fastica(sig_epoch');
    
    i_icasig{i} = icasig_temp;
    i_mixing_mat{i} = mixing_mat_temp;
    i_sep_mat{i} = sep_mat_temp;
    
end


% visualize the trial by trial ICA components

% 3 for 25 ms time after, 2 for 100 ms after

num_modes_subtract = 2;

numInt = min(size(icasig_temp,1),10);

for j = 1:size(dataIntTime,3)
    figure
    for i = 1:numInt
        subplot(numInt,1,i)
        plot(i_icasig{j}(i,:))
        title(['ICA component # , ', num2str(i)])
    end
    subtitle(['Trial # ', num2str(j)])
    
end

% subtract each one of these components

subtracted_sig_cell = {};

for i = 1:size(dataIntTime,3)
    
    
    combined_ica_recon = (i_mixing_mat{i}(:,1:num_modes_subtract)*i_icasig{i}(1:num_modes_subtract,:))';
    
    subtracted_sig_ICA_temp = dataIntTime(:,:,i) - combined_ica_recon;
    
    subtracted_sig_cell{i} = subtracted_sig_ICA_temp;
    subtracted_sig_matrix(:,:,i) = subtracted_sig_ICA_temp;
    
    figure
    plot(subtracted_sig_ICA_temp(:,channelInt))
    hold on
    plot(dataIntTime(:,channelInt,i))
    title(['Number of ICA modes subtracted = ', num2str(num_modes_subtract)])
    legend({'subtracted signal','original signal'})
    title(['Channel ', num2str(channelInt), ' Trial ', num2str(i) ])
    
    figure
    plot(subtracted_sig_ICA_temp(:,channelInt))
    title(['Subtracted Signal for ', num2str(num_modes_subtract), ' ICA modes, Channel ', num2str(channelInt), ' Trial ', num2str(i)])
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% do it again!

i_icasig = {};
i_mixing_mat = {};
i_sep_mat = {};
scale_fac = 100;

num_modes_subtract = 2;

for i = 1:size(dataIntTime,3)
    sig_epoch = scale_fac.*squeeze(subtracted_sig_matrix(:,:,i));
    [icasig_temp,mixing_mat_temp,sep_mat_temp] = fastica(sig_epoch');
    
    i_icasig{i} = icasig_temp;
    i_mixing_mat{i} = mixing_mat_temp;
    i_sep_mat{i} = sep_mat_temp;
    
end

numInt = min(size(icasig_temp,1),10);

for j = 1:size(dataIntTime,3)
    figure
    for i = 1:numInt
        subplot(numInt,1,i)
        plot(i_icasig{j}(i,:))
        title(['ICA component # , ', num2str(i)])
    end
    subtitle(['Trial # ', num2str(j)])
    
end

subtracted_sig_cell = {};

for i = 1:size(dataIntTime,3)
    
    
    combined_ica_recon = (i_mixing_mat{i}(:,1:num_modes_subtract)*i_icasig{i}(1:num_modes_subtract,:))';
    
    subtracted_sig_ICA_temp = dataIntTime(:,:,i) - combined_ica_recon./scale_fac;
    
    subtracted_sig_cell{i} = subtracted_sig_ICA_temp;
    subtracted_sig_matrix(:,:,i) = subtracted_sig_ICA_temp;
    
    figure
    plot(subtracted_sig_ICA_temp(:,channelInt))
    hold on
    plot(dataIntTime(:,channelInt,i))
    title(['Number of ICA modes subtracted = ', num2str(num_modes_subtract)])
    legend({'subtracted signal','original signal'})
    title(['Channel ', num2str(channelInt), ' Trial ', num2str(i) ])
    
    figure
    plot(subtracted_sig_ICA_temp(:,channelInt))
    title(['Subtracted Signal for ', num2str(num_modes_subtract), ' ICA modes, Channel ', num2str(channelInt), ' Trial ', num2str(i)])
    
end

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


numInt = min(size(icasig_temp,1),10);

for j = 1:size(dataIntTime,3)
    figure
    for i = 1:numInt
        subplot(numInt,1,i)
        plot(i_icasigS{j}(i,:))
        title(['ICA component # , ', num2str(i)])
    end
    subtitle(['Trial # ', num2str(j)])
    
end

%% subtract each one of these components

subtracted_sig_cellS_I = {};

% try 8 for 1000, 6 for 100 scale factor, 4 for 10 scale factor
num_modes_subtract = 8;
subtracted_sig_matrixS_I = [];

for i = 1:size(dataIntTime,3)
    
    
    combined_ica_recon = (i_mixing_matS{i}(:,1:num_modes_subtract)*i_icasigS{i}(1:num_modes_subtract,:))';
    
   % subtracted_sig_ICA_temp = dataIntTime(:,:,i) - combined_ica_recon./scale_factor;
        subtracted_sig_ICA_temp = dataIntTime(:,:,i) - combined_ica_recon./scale_factor;

    subtracted_sig_cellS_I{i} = subtracted_sig_ICA_temp;
    subtracted_sig_matrixS_I(:,:,i) = subtracted_sig_ICA_temp;
    
    figure
    plot(subtracted_sig_ICA_temp(:,channelInt))
    hold on
    plot(dataIntTime(:,channelInt,i))
    title(['Channel ', num2str(channelInt), ' Trial ', num2str(i), 'Number of ICA modes subtracted = ', num2str(num_modes_subtract)])
    legend({'subtracted signal','original signal'})
    
    figure
    plot(subtracted_sig_ICA_temp(:,channelInt))
    title(['Subtracted Signal for ', num2str(num_modes_subtract), ' ICA modes, Channel ', num2str(channelInt), ' Trial ', num2str(i)])
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% try it channel by channel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dataTrial = permute(dataIntTime,[1,3,2]);

i_icasigC = {};
i_mixing_matC = {};
i_sep_matC = {};

scale_factor = 10;


for i = 1:size(dataTrial,3)
    sig_matrix_channel = [];
    
    for j = 1:size(dataTrial,2)
        sig_epoch = scale_factor.*squeeze(dataTrial(:,j,i));
        
        sig_matrix_channel(j,:) = sig_epoch;
    end
    [icasig_temp,mixing_mat_temp,sep_mat_temp] = fastica(sig_matrix_channel);
    
    i_icasigC{i} = icasig_temp;
    i_mixing_matC{i} = mixing_mat_temp;
    i_sep_matC{i} = sep_mat_temp;
    
end


% visualize the trial by trial ICA components for a single channel

% channelInt = 58;

numInt = size(i_icasigC{channelInt},1);


figure
for j = 1:numInt
    subplot(numInt,1,j)
    plot(i_icasigC{channelInt}(j,:))
    title(['ICA component # , ', num2str(j)])
end
%subtitle(['Channel ' num2str(channelInt)])


%%
%subtract each one of these components

subtracted_sig_cellC_I = {};

% try 6 for 100 scale factor, 4 for 10 scale factor
num_modes_subtract = 1;

for i = 1:size(dataTrial,3)
    sig_matrix_channel = [];
    
    combined_ica_recon = (i_mixing_matC{i}(:,1:num_modes_subtract)*i_icasigC{i}(1:num_modes_subtract,:))';
    
    subtracted_sig_ICA_temp = dataTrial(:,:,i) - combined_ica_recon./scale_factor;
    
    subtracted_sig_cellC_I{i} = subtracted_sig_ICA_temp;
    subtracted_sig_matrixC_I(:,:,i) = subtracted_sig_ICA_temp;
%     
%     if i == channelInt
%         figure
%         
%         for j = 1:size(subtracted_sig_ICA_temp,2)
%             subplot(size(subtracted_sig_ICA_temp,2),1,j)
%             plot(subtracted_sig_matrixC_I(:,j,channelInt))
%             hold on
%             plot(dataTrial(:,j,channelInt))
%             title(['Trial ', num2str(j)])
%             
%         end
%         legend({'subtracted signal','original signal'})
%     end
%     %subtitle(['Channel ', num2str(channelInt), 'Number of ICA modes subtracted = ', num2str(num_modes_subtract)])
%     
    
    if i == channelInt
        figure
        
        for j = 1:size(subtracted_sig_ICA_temp,2)
            figure
            plot(subtracted_sig_matrixC_I(:,j,channelInt))
            hold on
            plot(dataTrial(:,j,channelInt))
            title(['Trial ', num2str(j)])
            
        end
        legend({'subtracted signal','original signal'})
    end
    %subtitle(['Channel ', num2str(channelInt), 'Number of ICA modes subtracted = ', num2str(num_modes_subtract)])
    
    
    
end

%%
% SVD part -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% good signal =     subtracted_sig_matrixS_I;
% notch the signal

% good channels? 


data_permuted  = permute(subtracted_sig_matrixS_I,[1,3,2]);

% stack the data

data_stacked = reshape(data_permuted,[size(data_permuted,1)*size(data_permuted,2),size(data_permuted,3)]);

ave_data = squeeze(mean(data_permuted,2));
figure
plot(ave_data(:,channelInt))

notched = notch(ave_data,[60 120 240 360],fs_data);

figure
plot(ave_data(:,channelInt))
hold on
plot(notched(:,channelInt))
legend({'Unfiltered','Filtered'});

%[t_post,fw,nCpost] = timeFrequencyAnalWavelet(
%sigCCEPs = ica_sigCCEPextract[data_stacked


%%

if(exist('ignore','var'))
    % ignore badChans for SVD
    goods = ones(size(data,2),1);
    goods(ignore) = 0;
    goods = logical(goods);
    
    %fullData(:,~goods) = 0;
    dataTrim = data(:,goods);
end

if(exist('goodChans','var'))
    goods = zeros(size(data,2),1);
    goods(goodChans) = 1;
    goods = logical(goods);
    dataTrim = data(:,goods);
end


% transpose data
dataSVD = dataTrim';
% data needs to be in m x n form, where m is the number of channels, and n
% is the time points (rows = sensors, columns = samples)

[u,s,v] = svd(dataSVD,'econ');

% SVD
[u,s,v] = svd(subtracted_sig_ICA','econ');

% look at temporal part - columns of v - indivividually and together
figure
% all together

modes = [1:3];
% get default colormap and use this to plot same as above
subDiag = diag(s);
subDiag = subDiag(1:length(modes));

% scale them by singular values for the combined plot, NOT for individual
numSubs = length(modes)+1;
subplot(numSubs,1,1)
plot(v(:,modes).*repmat(subDiag,[1 size(v,1)])')
title({'Temporal portion of the 3 modes', 'scaled by singular value'}), legend('show')

leg=cell(length(modes),1);
for i=1:length(modes)
    leg{i}=['mode ', num2str(i)];
end
legend(leg);

co = get(gca,'ColorOrder');

for i=1:length(modes)
    subplot(numSubs,1,i+1)
    plot(v(:,modes(i)), 'color',co(i,:));
    title(['Temporal portion of mode #: ', num2str(modes(i))])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimal Threshold 
m = size(dataSVD,1);
n = size(dataSVD,2);

% this y is just the singular values sorted
y = svd(dataSVD,'econ');

% this is the singular value cutoff 
coeffs = optimal_SVHT_coef(m/n,0);


% denoised matrix 

dataSVD_denoise = diag(dataSVD);
dataSVD_denoise ( dataSVD_denoise  < (optimal_SVHT_coef(m/n,0) * median(dataSVD_denoise )) ) = 0;
Xhat = u * diag(dataSVD_denoise ) * v';



