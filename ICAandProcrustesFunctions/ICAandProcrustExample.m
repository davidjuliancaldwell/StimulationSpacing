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
plotIt = false;
channelInt = 61;
bads = [12 52];
pre = -400;
post = 455;

[icaMat, icaCell,artMat,artCell,tNew] =  ica_artifact_remove(t,data,bads,pre,post,fs_data,scale_factor,numComponentsSearch,plotIt,channelInt);

%% procrustes 
plotIt = true;

pre = -5;
post = 40;
[d_mat,z_mat] = procrustes_metric(t,tNew,data,artMat,pre,post,bads,plotIt);
