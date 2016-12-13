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

figure
plot(data_stacked(:,58))


%% ICA THAT !

[icasig,mixing_mat,sep_mat] = fastica(data_stacked');

%%
figure
subplot(2,1,1)
plot(icasig(1,:))
subplot(2,1,2)
plot(icasig(2,:))

% reconstruct first ICA
first_ica_recon = (mixing_mat(:,1)*icasig(1,:))';
figure
plot(first_ica_recon(:,58))
%%
% subtract this from signal 

subtracted_sig = data_stacked - first_ica_recon;
figure
plot(subtracted_sig(:,58))
hold on
plot(data_stacked(:,58))
legend({'subtracted signal','original signal'})


%% try double subtraction 
combined_ica_recon = (mixing_mat*icasig)';

subtracted_sig_2ICA = data_stacked - combined_ica_recon;

figure
plot(subtracted_sig_2ICA(:,58))
hold on
plot(data_stacked(:,58))
legend({'subtracted signal','original signal'})
figure
plot(subtracted_sig_2ICA(:,58))

%% ICA again

[icasig2,mixing_mat2,sep_mat2] = fastica(subtracted_sig_2ICA');


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
scale_fac = 1000;
[icasig,mixing_mat,sep_mat] = fastica(scale_fac.*data_stacked');

%%

numInt = 10;

figure
for i = 1:numInt
    subplot(numInt,1,i)
    plot(icasig(i,:))
end

%%
% reconstruct first ICA
first_ica_recon = (mixing_mat(:,1)*icasig(1,:))';
figure
plot(first_ica_recon(:,58))
%%
% subtract this from signal 

subtracted_sig = data_stacked - first_ica_recon./scale_fac;
figure
plot(subtracted_sig(:,58))
hold on
plot(data_stacked(:,58))
legend({'subtracted signal','original signal'})


%% try double subtraction 
combined_ica_recon = (mixing_mat*icasig)';

subtracted_sig_2ICA = data_stacked - combined_ica_recon;

figure
plot(subtracted_sig_2ICA(:,58))
hold on
plot(data_stacked(:,58))
legend({'subtracted signal','original signal'})
figure
plot(subtracted_sig_2ICA(:,58))

%% ICA again

[icasig2,mixing_mat2,sep_mat2] = fastica(subtracted_sig_2ICA');
