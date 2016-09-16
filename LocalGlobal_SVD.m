%% Project local modes onto global modes
% So what we already did isn't in the correct format for the global mode,
% so first build the local X_k's and the global X

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
X_k = zeros(size(dataEpochedHigh,1), sum(goods), size(dataEpochedHigh, 3));
for i=1:size(dataEpochedHigh, 3)
    X_k(:,:,i) = dataEpochedHigh(:, goods, i);
end

X = reshape(X_k, [size(X_k,1), size(X_k,2)*size(X_k,3)]);

% SVD of the local matrices


% SVD of the global matrix



% etither give it SVDanalysis(.....,[],goodChans);
%or SVDanalysis(.......,badChans,[]);

fullData = true;
%[u,s,v] = SVDanalysis(dataStackedGood,stim_chans,fullData,badChans,[]);
modes = 1:3;
[u,s,v, dataSVD, goods] = SVDanalysis(dataStackedGood,stim_chans,fullData,[],goodChans, modes);

