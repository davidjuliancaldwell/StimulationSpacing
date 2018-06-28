%% Recursive rPCA TEST
% This runs tests to see if it's ok to shuffle the data and then run the
% rPCA
close all; clear all

% Add paths:
addpath('C:\Users\jcronin.CSENETID\Code\Matlab\Experiments\stimulation_spacing\resistivity_extraction')
addpath('C:\Users\jcronin.CSENETID\Code\Matlab\Experiments\stimulation_spacing\inexact_alm_rpca')
addpath('C:\Users\jcronin.CSENETID\Code\Matlab\Experiments\stimulation_spacing\inexact_alm_rpca\PROPACK')

% Load epoched data:
load('C:\Users\jcronin.CSENETID\Data\Subjects\3f2113\data\d6\Matlab\StimulationSpacing\1sBefore1safter\stim_12_52.mat')
L = size(dataEpoched, 1); % the number of samples in each epoch (will use this same amount for non-epoched data)

% Load original data (not epoched):
load('C:\Users\jcronin.CSENETID\Data\Subjects\3f2113\data\d6\Matlab\StimulationSpacing\3f2113_RawData\stim_12_52_raw.mat')

% Epoch the original data between stims
data = Wave.data;
data_fs = Wave.info.SamplingRateHz;
stim = Sing.data(:,1);
stim_fs = Sing.info.SamplingRateHz;

%% During development: run this and then can look at the stimulationExtractionBetween function (below); remove when done
fs_data = data_fs;
stimTrig = stim;
fs_stim = stim_fs;

%%
[epochedSignal] = stimulationExtractionBetween(data, data_fs, stim, stim_fs, L);


dataBetween = A_stimBetweenExtract;
dataEpoched = squeeze(getEpochSignal(data,sts-presamps,sts+postsamps+1));


% load('/Users/jcronin/Desktop/Data/stim_12_52.mat')
% addpath('/Users/jcronin/Desktop/Code/Experiments/stimulation_spacing/inexact_alm_rpca')
% addpath('/Users/jcronin/Desktop/Code/Experiments/stimulation_spacing/inexact_alm_rpca/PROPACK')

% define stimulation channels
stimChan1 = stim_chans(1);
stimChan2= stim_chans(2);
goodChans = false(1,size(dataEpochedHigh,2));
goodChans(1:64)= true;
goodChans(stim_chans) = false;
t_orig = t;

data(:,:,1) = squeeze(dataEpochedHigh(:, goodChans, 1))';
% dataUnshifted = data;

global_lambda = false;
gL = 0.11;
%% Recursively run the robust PCA:
iters = 1; % number of recursive iterations - DETERMINE THRU DATA!!
rng(12345) % set random number seed, so that this is repeatable 
thresh = 0.0001; % threshold for optimization of lambda through procrustes values
plotIt = false;

for count=1:iters % This will loop through the number of times that you want to recursively run the rPCA
    
    % Shifting/offsetting the data for each recursive iteration:
    randShift = randi(size(data,2), [1 size(data,1)]);
    dataShifted = zeros(size(data(:,:,1)));
    for i=1:size(data,1)
        temp = data(i,:,count);
        dataShifted(i,:) = temp([randShift(i):end 1:randShift(i)-1]);
    end
    
    % Optimize lambda value
    ii = 1;
    lambda(1) = 0.01;
    lambda(2) = 0.9;
    changeProcrustes = 1;
    
    while changeProcrustes>thresh
        
        % Apply robust PCA
        [R1, R2] = singleRPCA(dataShifted, lambda(ii), -1, -1, false);
        
        % Unshift the low-rank matrix, for procrustes analysis
        % TO DO - this unshifting may not be necessary - check (but timing
        % would be messed up...)
        LR = R1.';
        Sp = R2.';
        LR_Reshifted = zeros(size(LR));
        Sp_Reshifted = zeros(size(LR));
        for i=1:size(data,1)
            LR_Reshifted(i,[randShift(i):end 1:randShift(i)-1]) = LR(i,:);
            Sp_Reshifted(i,[randShift(i):end 1:randShift(i)-1]) = Sp(i,:);
        end
        lowRankTemp = LR_Reshifted';
%         dataS(:,:,count)=dataShifted;
        sparseTemp = Sp_Reshifted';
        
        % Procrustes analysis to update lambda
        % Procrustes on the stim artifact (compare data to sparse matrix)
        pre = -5;
        post = 40; 
        % TO DO: may need to change 'data' each time (should I really be
        % comparing this to the original data, or just the previous low
        % rank each time?)
        [d_mat,~] = procrustes_metric(t_orig, t_orig, data(:,:,count)', sparseTemp, pre, post, [], plotIt);
        mean_d_mat(ii,1) = mean(d_mat);

        
        figure
        subplot(2,1,1)
        data_orig = data(60,(t_orig>=pre & t_orig<=post),count);
        sparse_orig = sparseTemp(t_orig>=pre & t_orig<=post, 60, count);
        t_short = t_orig(t_orig>=pre & t_orig<=post);
        plot(t_short,data_orig), hold on, plot(t_short, sparse_orig')
        title(['Artifact, lambda = ', num2str(lambda(ii))])
        xlabel('time (ms)')

        
        
        % Procrustes on the data post stim (compare data to low rank matrix)
        pre = 5;
        post = 50; 
        [d_mat,~] = procrustes_metric(t_orig, t_orig, data(:,:,count)', lowRankTemp, pre, post, [], plotIt);
        mean_d_mat(ii,2) = mean(d_mat);
        
        subplot(2,1,2)
        data_orig = data(60,(t_orig>=pre & t_orig<=post),count);
        lr_orig = lowRankTemp(t_orig>=pre & t_orig<=post);
        t_short = t_orig(t_orig>=pre & t_orig<=post);
        plot(t_short, data_orig), hold on, plot(t_short, lr_orig')
        title(['post stim, lambda = ', num2str(lambda(ii))])
        xlabel('time (ms)')

        
        if ii>=2 % Now need to start defining new lambda values
            if ii==2
                lambda(ii+1) = (lambda(1)+lambda(2))/2;
            else
                temp = mean(mean_d_mat,2);
                [~, min_ind1] = min(temp);
%                 temp = mean_d_mat;
                temp(min_ind1) = 1;
                [~, min_ind2] = min(temp);
                lambda(ii+1) = (lambda(min_ind1)+lambda(min_ind2))/2;
            end
        end
        
        if ii>=2
            changeProcrustes = abs(mean_d_mat(ii)-mean_d_mat(ii-1));
        end
        
        ii = ii + 1;
        
    end
    
    % Plot the changes to the mean procrustes values and lambda
    figure
    scatter(lambda(1:end-1), mean_d_mat(:,1), 'b'), hold on
    scatter(lambda(1:end-1), mean_d_mat(:,2), 'r')
    scatter(lambda(1:end-1), mean(mean_d_mat,2), 'g')
    legend('Procrustes on sparse', 'Procrustes on low rank', 'Mean across')
    xlabel('lambda values'), ylabel('mean procrustes value across channels')
    title('Lambda vs. mean procrustes across 62 non-stim channels')
    
    % Now take the minimum mean procrustes d value and use that lambda for
    % this iteration.
    % UPDATE THIS: for now need to re-run the rPCA with that lambda
    [~, min_ind] = min(mean(mean_d_mat,2));
    lambdaFinal(count) = lambda(min_ind);
    
    % Run rPCA
    if global_lambda == true
        L = gL;
    else
        L = lambdaFinal(count);
        [R1, R2] = singleRPCA(dataShifted, L, -1, -1, false);
    end
    
    % Unshift the low-rank matrix, which will then be used recursively for
    % another rPCA
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

%% Procrustes with function
figure
subplot(2,1,1)
plot(t_orig, data(1,:,1)'), hold on, plot(t_orig, sparse(1,:,1))
subplot(2,1,2)
plot(t_orig, data(1,:,2)')

plotIt = true;
pre = 0;
post = 40;


data_orig = data(:,(t_orig>=pre & t_orig<=post),1);
sparse_orig = sparse(:,(t_orig>=pre & t_orig<=post),1);
figure
plot(data_orig(1,:)), hold on, plot(sparse_orig(1,:))




[d_mat,z_mat] = procrustes_metric(t_orig, t_orig, data(:,:,1)', sparse(:,:,1)', pre, post, [], plotIt);
disp('Finished Procrustes analysis')

%% Plot data
figure
subplot(2,2,1)
waterfall(CH, T, abs(dataPlot(:,:,1))), axis tight
xlabel('ch'), ylabel('time (ms)'), zlabel('ecog')
title('Original data')

subplot(2,2,2)
waterfall(CH, T, abs(dataPlot(:,:,2))), axis tight
xlabel('ch'), ylabel('time (ms)'), zlabel('ecog')
title(['First iter, \lambda = ', num2str(lambdaFinal(1))])

subplot(2,2,4)
waterfall(CH, T, abs(dataPlot(:,:,ceil(iters/2)+1))), axis tight
xlabel('ch'), ylabel('time (ms)'), zlabel('ecog')
title([num2str(ceil(iters/2)), 'th iter, \lambda = ', num2str(lambdaFinal(ceil(iters/2)-1))])

subplot(2,2,3)
waterfall(CH, T, abs(dataPlot(:,:,end))), axis tight
xlabel('ch'), ylabel('time (ms)'), zlabel('ecog')
title(['Last iter = #', num2str(iters), ', \lambda = ', num2str(lambdaFinal(end))])

%% Plot first and last low-rank and sparse matrices
figure
subplot(2,2,1)
waterfall(CH, T, abs(dataPlot(:,:,2))), axis tight
xlabel('ch'), ylabel('time (ms)'), zlabel('ecog')
title(['First iter Low-rank, \lambda = ', num2str(lambdaFinal(1))])

subplot(2,2,2)
waterfall(CH, T, abs(sparsePlot(:,:,1))), axis tight
xlabel('ch'), ylabel('time (ms)'), zlabel('ecog')
title(['First iter Sparse, \lambda = ', num2str(lambdaFinal(1))])

subplot(2,2,3)
waterfall(CH, T, abs(dataPlot(:,:,end))), axis tight
xlabel('ch'), ylabel('time (ms)'), zlabel('ecog')
title(['Last iter #', num2str(iters), ' Low-rank, \lambda = ', num2str(lambdaFinal(end))])

subplot(2,2,4)
waterfall(CH, T, abs(sparsePlot(:,:,end))), axis tight
xlabel('ch'), ylabel('time (ms)'), zlabel('ecog')
title(['Last iter #', num2str(iters), ' Sparse, \lambda = ', num2str(lambdaFinal(end))])

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
    title(['Ch.', num2str(chan_interest(i))])
    ylim([-1e-4, 1e-4])
end

subplot(length(chan),1,1)
title(['Original and low-rank data sets with: \lambda = ', num2str(lambdaFinal(1)),...
        ', iter 1; and \lambda = ', num2str(lambdaFinal(end)), ', iter ', num2str(iters), ...
        sprintf('\nCh. '), num2str(chan_interest(1))])
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
title(['Original and low-rank data sets, Ch. ', num2str(chan_interest), ' with: \lambda = ', num2str(lambdaFinal), ', iters = ', num2str(iters)])
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