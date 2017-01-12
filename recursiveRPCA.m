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
iters = 1000; % number of recursive iterations - DETERMINE THRU DATA!!
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
    lambda = 0.04;
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

%% Plot data
ch=1:62;
t=(0:length(data)-1)/fs_data*1000;
[T,CH]=meshgrid(t,ch);

figure
subplot(2,2,1)
waterfall(CH, T, abs(data(:,:,1))), axis tight
xlabel('ch'), ylabel('time (ms)'), zlabel('ecog')
title(['Original data, \lambda = ', num2str(lambda)])

subplot(2,2,2)
waterfall(CH, T, abs(data(:,:,2))), axis tight
xlabel('ch'), ylabel('time (ms)'), zlabel('ecog')
title('First iter')

subplot(2,2,4)
waterfall(CH, T, abs(data(:,:,ceil(iters/2)+1))), axis tight
xlabel('ch'), ylabel('time (ms)'), zlabel('ecog')
title([num2str(ceil(iters/2)), 'th iter'])

subplot(2,2,3)
waterfall(CH, T, abs(data(:,:,end))), axis tight
xlabel('ch'), ylabel('time (ms)'), zlabel('ecog')
title(['Last iter = #', num2str(iters)])

%% Plot first and last low-rank and sparse matrices
figure
subplot(2,2,1)
waterfall(CH, T, abs(data(:,:,2))), axis tight
xlabel('ch'), ylabel('time (ms)'), zlabel('ecog')
title(['First iter Low-rank, \lambda = ', num2str(lambda)])

subplot(2,2,2)
waterfall(CH, T, abs(sparse(:,:,1))), axis tight
xlabel('ch'), ylabel('time (ms)'), zlabel('ecog')
title(['First iter Sparse, \lambda = ', num2str(lambda)])

subplot(2,2,3)
waterfall(CH, T, abs(data(:,:,end))), axis tight
xlabel('ch'), ylabel('time (ms)'), zlabel('ecog')
title(['Last iter Low-rank, \lambda = ', num2str(lambda)])

subplot(2,2,4)
waterfall(CH, T, abs(sparse(:,:,end))), axis tight
xlabel('ch'), ylabel('time (ms)'), zlabel('ecog')
title(['Last iter Sparse, \lambda = ', num2str(lambda)])

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

%% Compare one channel
chan = 60;
figure
subplot(2,1,1)
plot(t, data(chan,:,1),'g')
hold on
plot(t, data(chan,:,2),'b')
plot(t, data(chan,:,end),'r')
legend('Original data', 'First low-rank matrix', 'Last low-rank matrix')
title(['Original and low-rank data sets, Ch. ', num2str(chan), ' with: \lambda = ', num2str(lambda), ', iters = ', num2str(iters)])
ylabel('ECoG (V)')

subplot(2,1,2)
plot(t, data(chan,:,2),'b')
hold on
plot(t, data(chan,:,end),'r')
legend('First low-rank matrix','Last low-rank matrix')
title(['First and last low-rank matrices, Ch. ', num2str(chan)])
ylabel('ECoG (V)')
xlabel('time (ms)')