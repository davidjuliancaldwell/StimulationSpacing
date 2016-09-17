%% Project local modes onto global modes
close all, clear all, clc

%% load a data file, must update below with the file you want to upload
% Or set allStimPairs to true;
fileName = 'stim_12_52';

allStimPairs = true;

% to add paths to all the subfolders
addpath(genpath(pwd))

%SPECIFIC ONLY TO DJC DESKTOP RIGHT NOW
% filePath = 'C:\Users\djcald\Google Drive\GRIDLabDavidShared\CSNE YSP2016\1sBefore1safter\'; 
% filePath = 'C:\Users\djcald\GoogleDrive\GRIDLabDavidShared\20f8a3\StimulationSpacingChunked\';

%SPECIFIC ONLY TO JAC DESKTOP RIGHT NOW
filePath = 'C:\Users\jcronin\Data\Subjects\3f2113\data\d6\Matlab\StimulationSpacing\1sBefore1safter\';

%SPECIFIC ONLY TO JAC Laptop RIGHT NOW
% filePath = '/Users/jcronin/Desktop/Data/3f2113/1sBefore1safter/';

if ~allStimPairs
    load([filePath, fileName, '.mat'])
end

%% Start building the data matrices to SVD
% post stim time window
post_begin = 5;
post_end = (450+post_begin);

% Build the local and global data matrices

% X_k = zeros(sum(goods), size(dataEpochedHigh,1), size(dataEpochedHigh, 3));
% for i=1:size(dataEpochedHigh, 3)
%     X_k(:,:,i) = dataEpochedHigh(:, goods, i)';
% end
% X = reshape(X_k, [size(X_k,1), size(X_k,2)*size(X_k,3)]);

filter_it = 'y';

if ~allStimPairs
    % Disregard bad channels
    prompt = {'what is the list of channels to IGNORE? e.g. 1:8,12 '};
    dlg_title = 'BadChannels';
    num_lines = 1;
    defaultans = {num2str(stim_chans)};
    answerChans = inputdlg(prompt,dlg_title,num_lines,defaultans);
    badChans = str2num(answerChans{1});

    goods = ones(size(dataEpochedHigh,2),1);
    goods(badChans) = 0;
    goods(65:end) = 0; % Only look at the grid for now
    goods = logical(goods);
    chansToStack = goods;
    X = dataStack(dataEpochedHigh,t,post_begin,post_end,chansToStack,[],[],[],fs_data,filter_it)';
    X_k = reshape(X, [size(X,1), size(X,2)/size(dataEpochedHigh, 3), size(dataEpochedHigh, 3)]);
    X_forSVD = X;
    X_k_forSVD = X_k;
    fileToLoad = 1; % Just need to know for later purposes that you only loaded in one file
else
    % Load in all datafiles one at a time and concatenate them to the X and
    % X_k matrices
    fileNamePrefix = 'stim_';
    stimPairs = {'4_60', '12_52', '20_44' '25_32' '26_31' '27_30' '28_29' '28_36'}; % There's probably a better way to do this, like reading in from the given path...
    fileToLoad = cellfun(@(s) strcat(fileNamePrefix, s), stimPairs, 'UniformOutput', false);
    X = [];
    allStimChans = [];
    for i=1:length(fileToLoad)
        load([filePath, fileToLoad{i}])
        chansToStack = 1:64;
        allStimChans = horzcat(allStimChans, stim_chans);
        X = horzcat(X, dataStack(dataEpochedHigh,t,post_begin,post_end,chansToStack,[],[],[],fs_data,filter_it)');
    end
    X_k = reshape(X, [size(X,1), size(X,2)/(size(dataEpochedHigh, 3)*length(fileToLoad)), size(dataEpochedHigh, 3)*length(fileToLoad)]);
    
    % Need to decide how do deal with the fact that we don't want to include
    % the stim channels when they were stimming, but also probably shouldn't
    % discard all 15 stim channels
    goods = ones(size(X,1),1);
    goods(unique(allStimChans)) = 0;
    goods = logical(goods);
    
%     X_forSVD = X(goods,:);
%     X_k_forSVD = X_k(goods,:,:);

    % Or try keeping all channels
    X_forSVD = X;
    X_k_forSVD = X_k;
    goods = ones(size(X,1),1);
end

%% SVD of the local matrices
% Transpose data: data needs to be in m x n form, where m is the number of
% channels, and n is the time points (rows = sensors, columns = samples)
uL=zeros(sum(goods),sum(goods),size(X_k,3));
sL=zeros(size(uL));
vL=zeros(size(X_k,2),sum(goods),size(X_k,3));

for i=1:size(X_k,3)
    [uL(:,:,i), sL(:,:,i), vL(:,:,i)] = svd(X_k_forSVD(:,:,i),'econ');
end

% SVD of the global matrix
[uG, sG, vG] = svd(X_forSVD, 'econ');
figure
plot(diag(sG),'ko','Linewidth',[2])

modes=1:3;
SVDplot(uG, sG, vG, false, [], modes)


%% Plot projections
modes=3:5;
a=zeros(1,length(modes));

figure
colorIncrement=0.1;
labels = cell(1, size(dataEpochedHigh,3));
for pairs=1:length(fileToLoad)
    gcolor=1.0; % this is to control the color of the line
    for ii=1:size(dataEpochedHigh,3)
        ind=1;
        for i=modes
            a(ind) = uG(:,i).'*uL(:,i,ii+(pairs-1)*size(dataEpochedHigh,3));
            ind=ind+1;
        end
        scatter3(a(1),a(2),a(3), [], [0.0 gcolor 1.0],'filled')
        hold on
        gcolor=gcolor-colorIncrement;
        labels{ii}=['Pulse ', num2str(ii)];
    end
end
xlabel(['Mode ' num2str(modes(1))])
ylabel(['Mode ' num2str(modes(2))])
zlabel(['Mode ' num2str(modes(3))])
legend(labels)

if allStimPairs
    title('All stim pairs: 8 pairs*10 stim pulses, local mode onto global')
else
    title([strrep(fileName, '_', '\_'), ': each point is one stim pulse local mode onto global mode'])
end

