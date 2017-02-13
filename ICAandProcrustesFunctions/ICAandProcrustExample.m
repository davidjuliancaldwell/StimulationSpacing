%% DJC 1-26-2017 - script demonstrating ica_artifact_remove.m and procrustes_metric.m
% this is a script demonstrating how to use the ica_artifact_remove.m and
% procrucrstes_metric.m scripts

%%
% clear work space
close all;clear all;clc

% change directory to whatever working directory is desired
cd 'C:\Users\djcald\SharedCode\stimulation_spacing'

% make sure FastICA_25 is on the path
%%
% load data
% DJC desktop

load('C:\Users\djcald\Google Drive\GRIDLabDavidShared\CSNE YSP 2016\1sBefore1safter\stim_12_52.mat')
%load('C:\Users\djcald\Google Drive\GRIDLabDavidShared\CSNE YSP 2016\1sBefore1safter\stim_28_29.mat')


%%

% run ICA

% some starting parameters
data = dataEpochedHigh(:,1:64,:);
scale_factor = 1000;
numComponentsSearch = 10;
plotIt = true;
channelInt = 60;
bads = [12 52];
pre = -400;
post = 455;

[icaMat, icaCell,artMat,artCell,tNew] =  ica_artifact_remove(t,data,bads,pre,post,fs_data,scale_factor,numComponentsSearch,plotIt,channelInt);

%% procrustes
plotIt = true;

pre = -5;
post = 40;
[d_mat,z_mat] = procrustes_metric(t,tNew,data,artMat,pre,post,bads,plotIt);

%% try time frequency analysis on the data?

signal = icaMat;
pre_begin = -400;
pre_end = -5;
post_begin = 5;
post_end = 455;
filter_it = 'y';
[ncPost_m,t_post,fw] = make_waveletMatrix(signal,tNew,pre_end,pre_begin,post_begin,post_end,fs_data,filter_it,[]);

%% plot example

figure
trialInt = 5;
nCPost = squeeze(ncPost_m(:,:,channelInt,trialInt));

% set colormap using cbrewer
CT = cbrewer('div','RdBu',11);
% flip it so red is increase, blue is down
CT = flipud(CT);
imagesc(t_post, fw, nCPost); axis xy
xlabel('Time in ms')
ylabel('Frequency (Hz)')
colorbar
title('Wavelet analysis for Post Stim - Normalized to Pre')
colormap(CT);
set_colormap_threshold(gcf, [-1 1], [-10 10], [.5 .5 .5])
set(gca,'fontsize',14)

figure
plot(tNew,squeeze(icaMat(:,channelInt,trialInt)))
set(gca,'fontsize',14)
title('Artifact Reduced Signal for One Trial')
xlabel('Time in ms')
ylabel(['Signal \muV'])
