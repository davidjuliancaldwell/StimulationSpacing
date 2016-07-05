%% 7/4/2016 - spectral analysis script DJC


%% This is to plot the time series and do the FFT of Pre and Post 

% example channel for stims 28_29, 21.
% set IDX to be whatever we want to look at channel wise. Can extend to
% whole matrix/data set later
idx = 21;
filter_it = input('notch filter? input "yes" or "no"','s');

% get the low signal
sig = mean(dataEpochedLow(:,idx,:),3);

% extract pre
t_pre = t(t<0);

% extract post 
t_post = t(t>5);

% filter it?
if strcmp(filter_it,'yes')
sig_pre = notch(sig(t<0),[60 120 180 240],fs_data);
sig_post = notch(sig(t>5),[60 120 180 240],fs_data);
else
    sig_pre = sig(t<0);
    sig_post = sig(t>5);

end
    

figure

% plot pre and post 
% plotting/etc is in spectralAnalysis function 
[f_pre,P1_pre] = spectralAnalysis(fs_data,t_pre,sig_pre); 

[f_post,P1_post] = spectralAnalysis(fs_data,t_post,sig_post);

% get middle signal 
sig = mean(dataEpochedMid(:,idx,:),3);

if strcmp(filter_it,'yes')
sig_post = notch(sig(t>5),[60 120 180 240],fs_data);
else
    sig_post = sig(t>5);

end

sig_post = sig(t>5);
sig_post = notch(sig(t>5),[60 120 180 240],fs_data);


[f_post,P1_post] = spectralAnalysis(fs_data,t_post,sig_post);

% get high signal 

sig = mean(dataEpochedHigh(:,idx,:),3);

if strcmp(filter_it,'yes')
sig_post = notch(sig(t>5),[60 120 180 240],fs_data);
else
    sig_post = sig(t>5);

end

[f_post,P1_post] = spectralAnalysis(fs_data,t_post,sig_post);

% set legend 
legend({'pre','low','mid','high'})

%% This is to stack the data so we can do SVD and DMD 
%close all;

% example channel for stim_28_29 
idx = 21;

% this selects all of the High data more than 5 ms after the stim (let's
% ignore stim for now)

dataNoStim = dataEpochedHigh(t>5,:,:);

% get an example Channel 
dataStacked = dataNoStim(:,idx,:);
% stack that example channel 
dataStacked = dataStacked(:);
figure
plot(dataStacked);

% reshift the data so we can do a vector operation to stack all of it like
% we did for that one example 
data_permuted  = permute(dataNoStim,[1,3,2]);

% stack the data

data_stacked = reshape(data_permuted,[size(data_permuted,1)*size(data_permuted,2),size(data_permuted,3)]);
figure
% compare the below plot to that example stacked channel above to confirm 
plot(data_stacked(:,idx))

% make a vector of all of the channels we have 
goods = ones(80,1);

% pick the ones to ignore
bads = [28,29,72:80];
goods(bads) = 0;

% make a logical matrix
goods = logical(goods);

% select the good channels
dataStackedGood = data_stacked(:,goods);

% decide if we want to filter it 
notch_stacked = input('notch the data? "yes" or "no"','s');
if strcmp(notch_stacked,'yes')
    dataStackedGood = notch(dataStackedGood,[60 120 180 240],fs_data);
    figure
    % plot the filtered data for a sanity check 
    plot(dataStackedGood(:,idx));
end

%% This is doing a SVD of our data matrix 
% looks at the first 3 modes in space, time 
[u,s,v] = svd(dataStackedGood','econ');

figure
plot(diag(s),'ko','Linewidth',[2])
% to get percentage in mode
subplot(2,1,1) % plot normal
plot(diag(s)/sum(diag(s)),'ko','Linewidth',[2])
title('singular values, fractions')
set(gca,'fontsize',14)

subplot(2,1,2) % plot semilog
semilogy(diag(s)/sum(diag(s)),'ko','Linewidth',[2])
title('singular values, fractions, semilog plot')
set(gca,'fontsize',14)

% look at the modes in space 
figure
x = [1:size(dataStackedGood,2)];
plot(x,u(:,1:3),'Linewidth',[2])
title('mode spatial locations'), legend('show')
legend({'mode 1','mode 2','mode 3'});


% look at temporal part - columns of v
figure

plot(v(:,1:3),'Linewidth',[2])
title('Temporal portion of the 3 modes'), legend('show')
legend({'mode 1','mode 2','mode 3'});

% BELOW THIS IS CURRENTLY NOT FUN
%% dmd - this is trying to do DMD - I don't think there's much useful from here until we talk to them

Xraw = dataStackedGood';
dt = 1/fs_data;

% added in dt optional argument, dt is our sampling frequency
% added in number of stacks. Using 5 for right now. The paper talks about
% hn > 2m, where h is the stack number, n is the number of channels, and m
% is time snapshots

% r sets rank truncation
[Phi, mu, lambda, diagS, x0] = DMD(Xraw,'dt',dt,'nstacks',50);

% look at modes from SVD augmented data matrix
% look at diagonal of matrix S - singular values

figure
plot(diagS(1:5),'ko','Linewidth',[2])

% to get percentage in mode
subplot(2,1,1) % plot normal
plot(diagS(1:5)/sum(diagS(1:5)),'ko','Linewidth',[2])
title('singular values, fractions')
set(gca,'fontsize',14)

subplot(2,1,2) % plot semilog
semilogy(diagS(1:5)/sum(diagS(1:5)),'ko','Linewidth',[2])
title('singular values, fractions, semilog plot')
set(gca,'fontsize',14)

%% more DMD reconstruction - similarly, not of much utility yet. 

t = length(Xraw');

% M is number of points to reconstrunct
M = t;

% add in time delay of 10 ms for plotting correctly
t2 = [1:M]/fs_data;

[Xhat z0] = DMD_recon(Phi, lambda, x0, M);

figure
x = [1:min(size(Xraw))];

% only plot 1:numChans of Xhat, there rest are redundant?
XhatUnique = Xhat([1:min(size(Xraw))],:);
surf(x,t2,1e6*real(XhatUnique).'), shading interp, colormap('parula')
title({'Reconstruction and Projection of Signal' 'with Dynamic Mode Decomposition'})
set(gca,'Fontsize',[14])
xlabel('Channel')
ylabel('Time s')
zlabel('Magnitude (\muV)');
%ylim(1e3*[0.01 max(t2)])


% optional argument plotit, either a 0 or a 1 , 1 plots the spectrum
[f P] = DMD_spectrum(Phi, mu, 'plotit',[1]);
set(gca,'fontsize',[14])
title('DMD Spectrum')

[f P] = DMD_spectrum(Phi, mu, 'plotit',[1]);
set(gca,'fontsize',[14])
title({'DMD Spectrum','Limited Frequency range'})
xlim([0 200])
% add in polarplot

if exist('polarplot')
    figure
    polarplot(mu,'*');
    title('Polar plot of \omega values from DMD')
    set(gca,'fontsize',[14])
end
