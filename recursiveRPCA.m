%% Recursive rPCA
clear all

% Load data:
load('C:\Users\jcronin\Data\Subjects\3f2113\data\d6\Matlab\StimulationSpacing\1sBefore1safter\stim_12_52.mat')
addpath('C:\Users\jcronin\Code\Matlab\Experiments\stimulation_spacing\inexact_alm_rpca')
addpath('C:\Users\jcronin\Code\Matlab\Experiments\stimulation_spacing\inexact_alm_rpca\PROPACK')

% define stimulation channels
stimChan1 = stim_chans(1);
stimChan2= stim_chans(2);
goodChans = false(1,size(dataEpochedHigh,2));
goodChans(1:64)= true;
goodChans(stim_chans) = false;

data(:,:,1) = squeeze(dataEpochedHigh(:, goodChans, 1))';
% dataUnshifted = data;


%% Recursively run the robust PCA:
iters = 1; % number of recursive iterations - DETERMINE THRU DATA!!
rng(12345) % set random number seed, so that this is repeatable 

for count=1:iters
    % Shifting/offsetting the data
    randShift = randi(size(data,2), [1 size(data,1)]);
    dataShifted = zeros(size(data(:,:,1)));
    for i=1:size(data,1)
        temp = data(i,:,count);
        dataShifted(i,:) = temp([randShift(i):end 1:randShift(i)-1]);
    end
    
    % Apply robust PCA
    lambda = 1;
    [R1, R2] = singleRPCA(dataShifted, lambda, -1, -1, false);
    
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
    data(:,:,count+1) = LR_Reshifted;
    dataS(:,:,count)=dataShifted;
    sparse(:,:,count) = Sp_Reshifted;
end
disp('rPCA loop finished')

% Create plot data matrices
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

%% Plot data
figure
subplot(2,2,1)
waterfall(CH, T, abs(dataPlot(:,:,1))), axis tight
xlabel('ch'), ylabel('time (ms)'), zlabel('ecog')
title(['Original data, \lambda = ', num2str(lambda)])

subplot(2,2,2)
waterfall(CH, T, abs(dataPlot(:,:,2))), axis tight
xlabel('ch'), ylabel('time (ms)'), zlabel('ecog')
title('First iter')

subplot(2,2,4)
waterfall(CH, T, abs(dataPlot(:,:,ceil(iters/2)+1))), axis tight
xlabel('ch'), ylabel('time (ms)'), zlabel('ecog')
title([num2str(ceil(iters/2)), 'th iter'])

subplot(2,2,3)
waterfall(CH, T, abs(dataPlot(:,:,end))), axis tight
xlabel('ch'), ylabel('time (ms)'), zlabel('ecog')
title(['Last iter = #', num2str(iters)])

%% Plot first and last low-rank and sparse matrices
figure
subplot(2,2,1)
waterfall(CH, T, abs(dataPlot(:,:,2))), axis tight
xlabel('ch'), ylabel('time (ms)'), zlabel('ecog')
title(['First iter Low-rank, \lambda = ', num2str(lambda)])

subplot(2,2,2)
waterfall(CH, T, abs(sparsePlot(:,:,1))), axis tight
xlabel('ch'), ylabel('time (ms)'), zlabel('ecog')
title(['First iter Sparse, \lambda = ', num2str(lambda)])

subplot(2,2,3)
waterfall(CH, T, abs(dataPlot(:,:,end))), axis tight
xlabel('ch'), ylabel('time (ms)'), zlabel('ecog')
title(['Last iter Low-rank, \lambda = ', num2str(lambda)])

subplot(2,2,4)
waterfall(CH, T, abs(sparsePlot(:,:,end))), axis tight
xlabel('ch'), ylabel('time (ms)'), zlabel('ecog')
title(['Last iter Sparse, \lambda = ', num2str(lambda)])

%% Check a few channels
chan_interest = [20 51 60 61]; % This is the true ECoG channel. You do not need to account for the stim channels or any bad channels
chan=zeros(size(chan_interest));
for i=1:length(chan)
    chan(i) = sum(goodChans(1:chan_interest(i))); % This will pull out the stim channels or any other channels that weren't plotted
end

figure
for i=1:length(chan)
    subplot(length(chan),1,i)
    plot(t, data(chan(i),:,1),'g', 'LineWidth', 1.2)
    hold on
    plot(t, data(chan(i),:,2),'b')
    plot(t, data(chan(i),:,end),'r')
    title(['Ch.', num2str(chan_interest(i)), ' with: \lambda = ', num2str(lambda), ', iters = ', num2str(iters)])
    ylim([-1e-4, 1e-4])
end

subplot(length(chan),1,1)
title(['Original and low-rank data sets, Ch. ', num2str(chan_interest(1)), ' with: \lambda = ', num2str(lambda), ', iters = ', num2str(iters)])
legend('Original data', 'First low-rank matrix', 'Last low-rank matrix')
ylabel('ECoG (V)')
%% Compare one channel
chan_interest = 60; % This is the true ECoG channel. You do not need to account for the stim channels or any bad channels
chan = sum(goodChans(1:chan_interest)); % This will pull out the stim channels or any other channels that weren't plotted
figure
subplot(2,1,1)
plot(t, data(chan,:,1),'g')
hold on
plot(t, data(chan,:,2),'b')
plot(t, data(chan,:,end),'r')
legend('Original data', 'First low-rank matrix', 'Last low-rank matrix')
title(['Original and low-rank data sets, Ch. ', num2str(chan_interest), ' with: \lambda = ', num2str(lambda), ', iters = ', num2str(iters)])
ylabel('ECoG (V)')

subplot(2,1,2)
plot(t, data(chan,:,2),'b')
hold on
plot(t, data(chan,:,end),'r')
legend('First low-rank matrix','Last low-rank matrix')
title(['First and last low-rank matrices, Ch. ', num2str(chan_interest)])
ylabel('ECoG (V)')
xlabel('time (ms)')

%% Procrustes values
d = zeros(1, size(data,1));
for i=1:size(data, 1)
    [d(i),Z,transform] = procrustes(data(i,:,1)',sparse(i,:,1)');
end

%% Procrutes, like David's
pre_procrust = -5;
post_procrust = 40;

data_orig =  dataIntTime((timeVec>pre_procrust & timeVec<post_procrust),:,:);

artifact = recon_artifact_matrix((timeVec>pre_procrust & timeVec<post_procrust),:,:);

procrustes_mat = zeros(size(artifact));


% channel by channel, trial by trial
for i = 1:numTrials
    
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

%% plot procrustes values
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
average_proc_trial = mean(d_mat,1)


figure
hold on 
for i = 1:numChans
    
    scatter(repmat(i,size(d_mat,2),1),d_mat(i,:))
    
end
xlabel(['Channel Number'])
ylabel(['Procrustes goodness of fit'])
title(['Procrustes Distance Metric for Artifact vs. Original Signal'])
sprintf('Average procrustes GOF for each channel')

average_proc_channel = mean(d_mat,2)

a = d_mat(:);





%% Plot the shifted data
% figure
% subplot(2,2,1)
% waterfall(CH, T, abs(dataS(:,:,1))), axis tight
% xlabel('ch'), ylabel('time (ms)'), zlabel('ecog')
% title('1. Shifted original data')
% 
% subplot(2,2,2)
% waterfall(CH, T, abs(dataS(:,:,2))), axis tight
% xlabel('ch'), ylabel('time (ms)'), zlabel('ecog')
% title('2. Shifted first iter')
% 
% subplot(2,2,3)
% waterfall(CH, T, abs(dataS(:,:,6))), axis tight
% xlabel('ch'), ylabel('time (ms)'), zlabel('ecog')
% title('3. Shifted 6th iter')
% 
% subplot(2,2,4)
% waterfall(CH, T, abs(dataS(:,:,10))), axis tight
% xlabel('ch'), ylabel('time (ms)'), zlabel('ecog')
% title('4. Shifted last iter')