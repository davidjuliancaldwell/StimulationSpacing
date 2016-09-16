%% Project local modes onto global modes
%% load a data file, must update below with the file you want to upload
% Or set allStimPairs to true;
fileName = 'stim_28_29';

allStimPairs = false;

% to add paths to all the subfolders
addpath(genpath(pwd))

%SPECIFIC ONLY TO DJC DESKTOP RIGHT NOW
% load(['C:\Users\djcald\Google Drive\GRIDLabDavidShared\CSNE YSP 2016\1sBefore1safter\', file Name, '.mat'])
% load(['C:\Users\djcald\Google Drive\GRIDLabDavidShared\20f8a3\StimulationSpacingChunked\', fileName, '.mat'])

%SPECIFIC ONLY TO JAC DESKTOP RIGHT NOW
load(['C:\Users\jcronin\Data\Subjects\3f2113\data\d6\Matlab\StimulationSpacing\1sBefore1safter\', fileName, '.mat'])

%SPECIFIC ONLY TO JAC Laptop RIGHT NOW
% load(['/Users/jcronin/Desktop/Data/3f2113/1sBefore1safter/', fileName, '.mat'])


%% Start building the data matrices to SVD
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

% Build the local and global data matrices
X_k = zeros(sum(goods), size(dataEpochedHigh,1), size(dataEpochedHigh, 3));
for i=1:size(dataEpochedHigh, 3)
    X_k(:,:,i) = dataEpochedHigh(:, goods, i)';
end

X = reshape(X_k, [size(X_k,1), size(X_k,2)*size(X_k,3)]);

%% SVD of the local matrices
% Transpose data: data needs to be in m x n form, where m is the number of
% channels, and n is the time points (rows = sensors, columns = samples)

uL=zeros(sum(goods),sum(goods),size(X_k,3));
sL=zeros(size(uL));
vL=zeros(size(X_k,2),sum(goods),size(X_k,3));

for i=1:size(X_k,3)
    dataSVD = X_k(:,:,i)';
    fullData = true;
    [uL(:,:,i), sL(:,:,i), vL(:,:,i)] = svd(X_k(:,:,i),'econ');
%     modes = 5:8;
%     SVDplot(uL(:,:,i), sL(:,:,i), vL(:,:,i), fullData, goods, modes);
end

% SVD of the global matrix
[uG, sG, vG] = svd(X, 'econ');

%% Plot projections
modes=4:6;
a=zeros(1,length(modes));

figure
for ii=1:size(X_k,3)
    ind=1;
    for i=modes
        a(ind) = uG(:,i).'*uL(:,i,ii);
        ind=ind+1;
    end
    scatter3(a(1),a(2),a(3), 'filled')
    hold on
end
xlabel(['Mode ' num2str(modes(1))])
ylabel(['Mode ' num2str(modes(2))])
zlabel(['Mode ' num2str(modes(3))])

if allStimPairs
else
    title([strrep(fileName, '_', '\_'), ': each point is one stim pulse local mode onto global mode'])
end

