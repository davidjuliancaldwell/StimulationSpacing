%% Denoise Data
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

data = squeeze(dataEpochedHigh(:, goodChans, 1))';
dataUnshifted = data;
% n=200;


%% Try shifting/offsetting the data
rng(12345)
randShift = randi(size(data,2), [1 size(data,1)]);
for i=1:size(data,1)
   temp = data(i,:);
   dataShifted(i,:) = temp([randShift(i):end 1:randShift(i)-1]);
end

data = dataShifted;

%% Plot data
ch=1:62;
t=(0:length(data)-1)/fs_data*1000;
[T,CH]=meshgrid(t,ch);
figure(1)
subplot(2,2,1)
plot(t, data(30,:)), title('Ch. 30'), axis tight
subplot(2,2,2)
waterfall(CH, T, abs(data))
title('original data shifted')
axis tight
xlabel('ch'), ylabel('time (ms)'), zlabel('ecog')

%% Robust PCA
[U,S,V]=svd(dataShifted, 'econ');
[U_unshifted,S_unshifted,V_unshifted]=svd(dataUnshifted, 'econ');

%% Plot the singular values
figure(2)
subplot(2,2,1)
plot(diag(S_unshifted),'ko','Linewidth',[2])
title('singular values on unshifted data')
subplot(2,2,2)
plot(diag(S),'ko','Linewidth',[2])
title('singular values on shifted data')

%% Plot the temporal portions of the first 4 modes
figure(3)
subplot(2,1,1)
plot(t,V(:,1:4))
title('data'), xlabel('x'), ylabel('V')
legend('mode 1', 'mode 2', 'mode 3', 'mode 4')

%% Apply robust PCA
% the corrupt data matrix is decomposed into real and imaginary parts
% before applying the inexact_alm_rpca algorithm
ur=real(data);
% ui=imag(data);

% The parameter lambda used in the inexact_alm_rpca algorithm can be tuned
% to best separate the sparse from low-rank matrices
lambda=0.2; 
% See below for examples of changing lambda

% low-rank (R1) and sparse (R2) portions
[R1r, R2r] = inexact_alm_rpca(real(ur.'),lambda);
% [R1i, R2i] = inexact_alm_rpca(real(ui.'),lambda);

% the real and imaginary portions are put back together here to SVD the
% low-rank portion:
% R1 = R1r + 1i*R1i;
% R2 = R2r + 1i*R2i;
R1 = R1r;
R2 = R2r;

[Ur,Sr,Vr]=svd(R1.');

%% Unshift data (only if you shifted it)
% low-rank matrix:
LR = R1.';
Sp = R2.';
LR_Reshifted = zeros(size(LR));
Sp_Reshifted = zeros(size(LR));
for i=1:size(data,1)
%    LR_Reshifted([(end-(randShift(i)-2))]) = temp;
   LR_Reshifted(i,[randShift(i):end 1:randShift(i)-1]) = LR(i,:);
   Sp_Reshifted(i,[randShift(i):end 1:randShift(i)-1]) = Sp(i,:);
end

%% Plot unshifted data
figure
subplot(2,2,1)
waterfall(CH, T, abs(dataUnshifted)), axis tight
xlabel('ch'), ylabel('time (ms)'), zlabel('ecog')
title('1. Original data')

subplot(2,2,2)
waterfall(CH, T, abs(dataShifted)), axis tight
xlabel('ch'), ylabel('time (ms)'), zlabel('ecog')
title('2. Shifted data')

subplot(2,2,4)
% waterfall(CH, T, abs(LR)), axis tight
% xlabel('ch'), ylabel('time (ms)'), zlabel('ecog')
% title('3. Low-rank shifted data')
waterfall(CH, T, abs(Sp_Reshifted)), axis tight
xlabel('ch'), ylabel('time (ms)'), zlabel('ecog')
title('3. Sparse re-shifted data')

subplot(2,2,3)
waterfall(CH, T, abs(LR_Reshifted)), axis tight
xlabel('ch'), ylabel('time (ms)'), zlabel('ecog')
title('4. Low-rank re-shifted data')

%% Do another SVD on the low-rank, re-shifted (unshifted) data
[Ur_unshifted,Sr_unshifted,Vr_unshifted]=svd(LR_Reshifted);

%% Add the subplots for the RPCA results
figure(1)
subplot(2,2,3), plot(t, R1(:,30)), title('Ch. 30'), axis tight
subplot(2,2,4), waterfall(CH,T,abs(R1)'), colormap([0 0 0]), title('RPCA: low-rank')
axis tight

figure(2)
subplot(2,2,4)
plot(diag(Sr),'ko','Linewidth',[2])
title('RPCA: low-rank, shifted data')
subplot(2,2,3)
plot(diag(Sr_unshifted),'ko','Linewidth',[2])
title('RPCA: low-rank, re-shifted data')

figure(3)
subplot(2,1,2)
plot(t,Vr(:,1:4))
title('RPCA: low-rank'), xlabel('x'), ylabel('V')

figure
subplot(2,1,1)
plot(t,Vr_unshifted(:,1:4))
title('Temporal modes, data un-shifted'), xlabel('x'), ylabel('V')
legend('mode 1', 'mode 2', 'mode 3', 'mode 4')
subplot(2,1,2)
plot(t,Vr_unshifted(:,1:4))
title('RPCA: low-rank, re-shifted'), xlabel('x'), ylabel('V')

%% Plot the low-rank, sparse, and original data matrices
figure
subplot(2,2,1)
waterfall(CH, T, abs(dataUnshifted)), axis tight
xlabel('ch'), ylabel('time (ms)'), zlabel('ecog')
title('1. Original data')

subplot(2,2,2)
waterfall(CH, T, abs(dataShifted)), axis tight
xlabel('ch'), ylabel('time (ms)'), zlabel('ecog')
title('2. Shifted data')

subplot(2,2,3)
waterfall(CH, T, abs(R1.')), axis tight
xlabel('ch'), ylabel('time (ms)'), zlabel('ecog')
title('3. Low-rank data')

subplot(2,2,4)
waterfall(CH, T, abs(R2.')), axis tight
xlabel('ch'), ylabel('time (ms)'), zlabel('ecog')
title('4. Sparse')

%% Plot the temporal portions of the other modes
modes = 9:12;
figure(4)
subplot(2,1,1)
plot(t,V(:,modes))
title(['original data - temporal modes: ', num2str(modes)]), xlabel('x'), ylabel('V')
for i=1:length(modes)
    labels{i}=['mode ', num2str(modes(i))];
end
legend(labels)
subplot(2,1,2)
plot(t,Vr(:,modes))
title('robust PCA temporal modes'), xlabel('x'), ylabel('V')
legend(labels)

%% Plot the original shifted data versus R1
figure(5)
chan = 60;
% subplot(2,1,1)
plot(data(chan,:))
hold on 
% subplot(2,1,2)
plot(R1(:,chan))
legend('original data', 'rPCA low-rank component')

%% Plot some channels to check
figure
for i=1:62
    plot(dataUnshifted(i,:))
    hold on 
    plot(LR_Reshifted(i,:))
    scatter(randShift(i), dataUnshifted(i,randShift(i)),'k')
    legend('original data', 'rPCA low-rank component', 'shift point')
    title(['Original data vs. low-rank data (un-shifted), Ch. ', num2str(i)])
    hold off
    pause(0.5)
end

%% Check the boundary where channels had been shifted
figure
chan = 35;
plot(dataUnshifted(chan,:))
hold on
plot(LR_Reshifted(chan,:))
scatter(randShift(chan), dataUnshifted(chan,randShift(chan)),'k')

%% Run SVD on the re-shifted
[Ur,Sr,Vr]=svd(R1.');

%% Low-rank approximations using the SVD of the corrupt data matrix vs. the RPCA
% low-rank approximation of corrupt data matrix
figure(4)
for j=2:5
    ff=U2(:,1:j)*S2(1:j,1:j)*V2(:,1:j).';
    subplot(2,2,j-1), waterfall(x,t,abs(ff)), colormap([0 0 0])
    title(['rank= ', num2str(j)])
    xlabel('x'), ylabel('t'), zlabel('|u|')
end
subplot(2,2,1), title(['low-rank approximation of corrupt data matrix;' ' rank= ' num2str(2)])

% approximation of low-rank matrix from robust PCA
figure(5)
for j=2:5
    ff=U3(:,1:j)*S3(1:j,1:j)*V3(:,1:j).';
    subplot(2,2,j-1), waterfall(x,t,abs(ff)), colormap([0 0 0])
    title(['rank= ', num2str(j)])
    xlabel('x'), ylabel('t'), zlabel('|u|')
end
subplot(2,2,1), title(['low-rank approximation  from robust PCA;' ' rank= ' num2str(2)])

%% Plots for unshifted data (i.e. original data thru RPCA)
ch=1:62;
t=(0:length(data)-1)/fs_data*1000;
[T,CH]=meshgrid(t,ch);

figure
subplot(3,1,1)
waterfall(CH, T, abs(dataUnshifted)), axis tight
xlabel('ch'), ylabel('time (ms)'), zlabel('ecog')
title('Original data')

subplot(3,1,2)
waterfall(CH, T, abs(R2.')), axis tight
xlabel('ch'), ylabel('time (ms)'), zlabel('ecog')
title('RPCA Sparse, from never shifted data')

subplot(3,1,3)
waterfall(CH, T, abs(R1.')), axis tight
xlabel('ch'), ylabel('time (ms)'), zlabel('ecog')
title('RPCA Low-rank, from never shifted data')

% Then plot a comparison plot for one chan:
figure
chan = 21;
plot(data(chan,:))
hold on 
plot(R1(:,chan))
legend('original data', 'rPCA low-rank component')


%% Old stuff below here:
% DBS data is organized as: time x channels x epochs
data = ECoG_sep{4};
epoch=1;

figure(1)
for i=1:8
    chan=i;
    subplot(4,2,i)
    plot(DBS_sep{4}(:,chan,epoch))
    title(['DBS electrode ', num2str(chan)]) 
    xlabel('samples'), ylabel('V')
end

figure(2)
for i=1:8
    chan=i;
    subplot(4,2,i)
    plot(data(:,chan,epoch))
    title(['ECoG electrode ', num2str(chan)]) 
    xlabel('samples'), ylabel('V')
end

figure(3)
for i=1:8
    chan=i+8;
    subplot(4,2,i)
    plot(data(:,chan,epoch))
    title(['ECoG electrode ', num2str(chan)]) 
    xlabel('samples'), ylabel('V')
end

%% DBS Spectrum
[amp, phase] = hilbAmp(data(:,:,1), [0 200], ECOG_fs);
figure
plot(amp)

figure
plot(abs(hilbert(data(:,:,1))))

figure
[pxx,f] = pwelch(data(:,:,1),ECOG_fs);
plot(pxx)


%% Mean subtract and visualize

%% Notch and visualize

dataF = highpass(data(:,:,epoch), 0.1, ECOG_fs, 2);
dataF = notch(dataF,[60 120 180 240],ECOG_fs, 2); % notched needs time x channels

% figure(2)
% hold on
figure
for i=1:8
    chan=i;
    subplot(4,2,i)
    plot(data(:,chan,epoch))
    hold on
    plot(dataF(:,chan),'r')
    title(['ECoG electrode ', num2str(chan)])
    xlabel('samples'), ylabel('V')
end
legend('original', 'highpassed and notched')

%% figure(3)
hold on
for i=1:8
    chan=i+8;
    subplot(4,2,i)
    plot(dataF(:,chan,epoch),'r')
    title(['ECoG electrode ', num2str(chan)]) 
    xlabel('samples'), ylabel('V')
end
legend('original', 'mean subtracted and notched')

