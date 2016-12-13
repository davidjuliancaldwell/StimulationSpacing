%% DJC -12-12-2016
% prototyping with ICA for stimulation extraction prior to data analysis
%%
close all;clear all;clc

%%

% DJC desktop

load('C:\Users\djcald\Google Drive\GRIDLabDavidShared\CSNE YSP 2016\1sBefore1safter\stim_12_52.mat')

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pre time window
pre_begin = -10;
pre_end = 0;
% post time window
post_begin = 5;
post_end = (100+post_begin);
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

[icasig,mixing_mat,sep_mat] = fastica(data_stacked');


%% subtract however many ICA components
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
title(['Number of ICA modes subtracted = ', num2str(num_modes_subtract)])
legend({'subtracted signal','original signal'})
subtitle(['Channel ' num2str(channelInt)])

figure
plot(subtracted_sig_ICA(:,58))
title(['Subtracted Signal for ', num2str(num_modes_subtract), ' ICA modes, Channel ', num2str(channelInt)])


%% SVD
[u,s,v] = svd(subtracted_sig_2ICA','econ');

% look at temporal part - columns of v - indivividually and together
figure
% all together

modes = [1:3]
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

%% try scaled ICA
scale_fac = 10;
[icasigS,mixing_matS,sep_matS] = fastica(scale_fac.*data_stacked');

num_modes_subtract = 6;

numInt = min(size(icasigS,1),10);

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
title(['Number of ICA modes subtracted = ', num2str(num_modes_subtract)])
legend({'subtracted signal','original signal'})
subtitle(['Channel ' num2str(channelInt)])

figure
plot(subtracted_sigS_ICA(:,58))
title(['Subtracted Signal for ', num2str(num_modes_subtract), ' ICA modes, Channel ', num2str(channelInt)])

%% Try single trials

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


%% visualize the trial by trial ICA components
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

%% subtract each one of these components

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
    subtitle(['Channel ', num2str(channelInt), ' Trial ', num2str(i) ])
    
    figure
    plot(subtracted_sig_ICA_temp(:,58))
    title(['Subtracted Signal for ', num2str(num_modes_subtract), ' ICA modes, Channel ', num2str(channelInt), ' Trial ', num2str(i)])
    
end


%% do it again! 

i_icasig = {};
i_mixing_mat = {};
i_sep_mat = {};
scale_fac = 100;

num_modes_subtract = 4 ;

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
    subtitle(['Channel ', num2str(channelInt), ' Trial ', num2str(i) ])
    
    figure
    plot(subtracted_sig_ICA_temp(:,58))
    title(['Subtracted Signal for ', num2str(num_modes_subtract), ' ICA modes, Channel ', num2str(channelInt), ' Trial ', num2str(i)])
    
end

%% single trial scaled 


i_icasigS = {};
i_mixing_matS = {};
i_sep_matS = {};

scale_factor = 10;


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

subtracted_sig_cell = {};

% try 6 for 100 scale factor, 4 for 10 scale factor 
num_modes_subtract = 4;

for i = 1:size(dataIntTime,3)
    
    
    combined_ica_recon = (i_mixing_matS{i}(:,1:num_modes_subtract)*i_icasigS{i}(1:num_modes_subtract,:))';
    
    subtracted_sig_ICA_temp = dataIntTime(:,:,i) - combined_ica_recon./scale_factor;
    
    subtracted_sig_cellS{i} = subtracted_sig_ICA_temp;
    subtracted_sig_matrixS(:,:,i) = subtracted_sig_ICA_temp;
    
    figure
    plot(subtracted_sig_ICA_temp(:,channelInt))
    hold on
    plot(dataIntTime(:,channelInt,i))
    title(['Number of ICA modes subtracted = ', num2str(num_modes_subtract)])
    legend({'subtracted signal','original signal'})
    subtitle(['Channel ', num2str(channelInt), ' Trial ', num2str(i) ])
    
    figure
    plot(subtracted_sig_ICA_temp(:,58))
    title(['Subtracted Signal for ', num2str(num_modes_subtract), ' ICA modes, Channel ', num2str(channelInt), ' Trial ', num2str(i)])
    
end



