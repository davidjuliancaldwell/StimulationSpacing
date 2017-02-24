close all; clear all

% Load data:
% load('C:\Users\jcronin\Data\Subjects\3f2113\data\d6\Matlab\StimulationSpacing\1sBefore1safter\stim_12_52.mat')
% addpath('C:\Users\jcronin\Code\Matlab\Experiments\stimulation_spacing\inexact_alm_rpca')
% addpath('C:\Users\jcronin\Code\Matlab\Experiments\stimulation_spacing\inexact_alm_rpca\PROPACK')

load('/Users/jcronin/Desktop/Data/stim_12_52.mat')
addpath('/Users/jcronin/Desktop/Code/Experiments/stimulation_spacing/inexact_alm_rpca')
addpath('/Users/jcronin/Desktop/Code/Experiments/stimulation_spacing/inexact_alm_rpca/PROPACK')

% define stimulation channels
stimChan1 = stim_chans(1);
stimChan2= stim_chans(2);
goodChans = false(1,size(dataEpochedHigh,2));
goodChans(1:64)= true;
goodChans(stim_chans) = false;
t_orig = t;

data(:,:,1) = squeeze(dataEpochedHigh(:, goodChans, 1))';
% dataUnshifted = data;

%% Run rPCA for multiple values of lambda:
lambda = [0.0001, 0.0001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
% iters = 5; % number of recursive iterations - DETERMINE THRU DATA!!
rng(12345) % set random number seed, so that this is repeatable 


% Shifting/offsetting the data
randShift = randi(size(data,2), [1 size(data,1)]);
dataShifted = zeros(size(data(:,:,1)));
for i=1:size(data,1)
    temp = data(i,:);
    dataShifted(i,:) = temp([randShift(i):end 1:randShift(i)-1]);
end


for count=1:length(lambda)
    % Apply robust PCA
%     lambda = 0.1;
    [R1, R2] = singleRPCA(dataShifted, lambda(count), -1, -1, false);
    
    % Unshift the low-rank matrix, which will then be used recursively for
    % another RPCA
    LR = R1.';
    Sp = R2.';
    LR_Reshifted = zeros(size(LR));
    Sp_Reshifted = zeros(size(LR));
    for i=1:size(data,1)
        LR_Reshifted(i,[randShift(i):end 1:randShift(i)-1]) = LR(i,:);
        Sp_Reshifted(i,[randShift(i):end 1:randShift(i)-1]) = Sp(i,:);
    end
    lowRank(:,:,count) = LR_Reshifted';
    dataS(:,:,count)=dataShifted';
    sparse(:,:,count) = Sp_Reshifted';
    
    

    
end
disp('rPCA loop finished')

%% Procrustes analysis to update lambda
plotIt = true;
pre = -5;
post = 40; 
[d_mat,z_mat] = procrustes_metric(t_orig, t_orig, lowRank, sparse, pre, post, [], plotIt);
mean_d_mat = mean(d_mat);

figure
scatter(lambda, mean_d_mat)

%% Create plot data matrices
full = 0; % Set this to 1 if you want to plot all of the data points, otherwise full=0 will just plot some data points to improve run time 

if full==1
    ch=1:62;
    t=(0:length(data)-1)/fs_data*1000;
    [T,CH]=meshgrid(t,ch);
    dataPlot=data;
    sparsePlot=sparse;
elseif full==0
    ch=1:62;
    step=5; % Plot data point every _step_ entries
    t=(0:step:length(data)-1)/fs_data*1000;
    [T,CH]=meshgrid(t,ch);
    dataPlot=data(:,1:step:end,:);    
    sparsePlot=sparse(:,1:step:end,:);
    t=(0:length(data)-1)/fs_data*1000;
end
