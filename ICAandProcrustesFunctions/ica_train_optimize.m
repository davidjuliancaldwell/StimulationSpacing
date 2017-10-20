function [evaluated] =  ica_train_optimize(scale_factor,data,stimChans,fs_data,meanSub,orderPoly)
%USAGE: function [subtracted_sig_matrixS_I, subtracted_sig_cellS_I] =  ica_artifact_remove(t,data,stimChans,pre,post,fs_data,scale_factor,numComponentsSearch,plotIt,channelInt)
%This function will perform the fast_ica algorithm upon a data set in the
%format of m x n x p, where m is samples, n is channels, and p is the
%individual trial. This is for trains of stimuli
%
% data = samples x channels x trials
% tTotal =  time vector
% stimChans = stimulation channels, or any channels to ignore
% pre = the time point at which to begin extracting the signal
% post = the time point at which to stop extracting the signal
% fs_data = sampling rate (Hz)
% scale_factor = scaling factor tp ensure the ICA algorithm functions
%       correctly
%numComponentsSearch = the number of ICA components to search through for
%       artifacts that meet a certain profile
% plotIt = plot it or not
% channelInt = plot a channel if interested
% REQUIRES FastICA algorithm in path

% set scale factor


% plot intermediate steps

if (~exist('meanSub','var'))
    meanSub = 0;
end

if (~exist('orderPoly','var'))
    orderPoly = 2;
end

freqFilter = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get stim channels, as we don't want to perform ICA on them

bads = [];
badTotal = [stimChans; bads];

% total channels
numChans = size(data,2);
% make logical good channels matrix to index
goods = zeros(numChans,1);

channelsOfInt = 1:numChans;

goods(channelsOfInt) = 1;
% set the goods matrix to be zero where bad channels are
goods(badTotal) = 0;
% make it logical
goods = logical(goods);

% make storage matrices
i_icasigS = {};
i_mixing_matS = {};
i_sep_matS = {};

% extract the data of interest
dataInt = data(:,goods,:);

% NOTE THIS IS DIFFERENT THAN BEFORE, WE WANT TO KEEP STIMULATION IN THERE
dataIntTime = dataInt;

if meanSub == 1
    for i = 1:size(dataIntTime,2)
        for j = 1:size(dataIntTime,3)
            data_int_temp = squeeze(dataIntTime(:,i,j));
            [p,s,mu] = polyfit((1:numel(data_int_temp))',data_int_temp,orderPoly);
            f_y = polyval(p,(1:numel(data_int_temp))',[],mu);
            
            % subtract poly fit
            dataIntTime(:,i,j) = data_int_temp - f_y;
            
            %dataIntTime = dataIntTime - repmat(mean(data,1),size(data,1),1);
            
        end
        %         figure;
        %         plot(f_y)
    end
end

numTrials = size(dataIntTime,3);

for i = 1:numTrials
    sig_epoch = scale_factor.*squeeze(dataIntTime(:,:,i));
    %[icasig_temp,mixing_mat_temp,sep_mat_temp] = fastica(sig_epoch','g','pow3','numOfIC',numComponentsSearch);
    [icasig_temp,mixing_mat_temp,sep_mat_temp] = fastica(sig_epoch','g','gauss','approach','symm','verbose','on');
    i_icasigS{i} = icasig_temp;
    i_mixing_matS{i} = mixing_mat_temp;
    i_sep_matS{i} = sep_mat_temp;
    
end


%% set ICA components that are like the artifact to zero (they occur near a certain time and have prominence)

% need to adjust this for case where it's close to zero but not quite
% equal?
numTrials = size(dataIntTime,3);

i_ica_art = {};
i_ica_mix_art = {};

% figure
% hold on

for i = 1:numTrials
    numICs = size(i_icasigS{i},1);
    start_index = 1;
    for j = 1:numICs
        
        % have to tune this
        [pk_temp_pos,locs_temp_pos] = findpeaks(i_icasigS{i}(j,:),fs_data,'MinPeakProminence',10);
        [pk_temp_neg,locs_temp_neg] = findpeaks(-1*i_icasigS{i}(j,:),fs_data,'MinPeakProminence',10);
        %
        
        %         findpeaks(-1*i_icasigS{i}(j,:),fs_data,'MinPeakProminence',10)
        %         findpeaks(i_icasigS{i}(j,:),fs_data,'MinPeakProminence',10)
        %         %
        
        % should be at least 10 peaks even at 100 Hz trains
        total_peaks = length(pk_temp_pos)+length(pk_temp_neg);
        
        [f,P1] = spectralAnalysisComp(fs_data,i_icasigS{i}(j,:));
        [maxi,ind] = max(P1(f>62));
        f_temp= f(f>62);
        rounded_f = round(f_temp(ind),-1);
        
        % DJC 4-17-2017 - add in 200 Hz stim peak
        if ((~isempty(pk_temp_pos) || ~isempty(pk_temp_neg)) && total_peaks > 10 && mod(rounded_f,200) == 0)
            
            i_ica_art{i}(start_index,:) = i_icasigS{i}(j,:);
            i_ica_mix_art{i}(:,start_index) = i_mixing_matS{i}(:,j);
            
            i_icasigS{i}(j,:) = 0;
            i_mixing_matS{i}(:,j) = 0; % according to ICA removes EEG
           % artifacts - just set rows of activation waveforms to be zero?
            
            start_index = start_index+1;
            
        else
            
            [f,P1] = spectralAnalysisComp(fs_data,i_icasigS{i}(j,:));
            [maxi,ind] = max(P1);
            rounded_f = round(f(ind),-1);
            
            %  200 Hz frequency content, 60 Hz frequency content (added in
            %  120, 180 4-11-2017
            %             if (mod(rounded_f,60) == 0 | mod(rounded_f,120) == 0 | mod(rounded_f,180) == 0)
            %                 %if mod(rounded_f,60) == 0 || mod(rounded_f,200) == 0
            %                 i_ica_kept{i}(start_index,:) = i_icasigS{i}(j,:);
            %                 i_ica_mix_kept{i}(:,start_index) = i_mixing_matS{i}(:,j);
            %                 start_index = start_index+1;
            %             end
            
        end
        %         if mod(rounded_f,200) == 0
        %             i_ica_kept{i}(start_index,:) = i_icasigS{i}(j,:);
        %             i_ica_mix_kept{i}(:,start_index) = i_mixing_matS{i}(:,j);
        %             start_index = start_index+1;
        %         end
        
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% this was the original
    %         % have to tune this
    %         [pk_temp_pos,locs_temp_pos] = findpeaks(i_icasigS{i}(j,:),fs_data,'MinPeakProminence',10);
    %         [pk_temp_neg,locs_temp_neg] = findpeaks(-1*i_icasigS{i}(j,:),fs_data,'MinPeakProminence',10);
    %         %
    %
    %         %         findpeaks(-1*i_icasigS{i}(j,:),fs_data,'MinPeakProminence',10)
    %         %         findpeaks(i_icasigS{i}(j,:),fs_data,'MinPeakProminence',10)
    %         %         %
    %
    %         % should be at least 10 peaks even at 100 Hz trains
    %         total_peaks = length(pk_temp_pos)+length(pk_temp_neg);
    %
    %         if ((~isempty(pk_temp_pos) || ~isempty(pk_temp_neg)) && total_peaks > 10)
    %
    %             i_ica_art{i}(start_index,:) = i_icasigS{i}(j,:);
    %             i_ica_mix_art{i}(:,start_index) = i_mixing_matS{i}(:,j);
    %
    %             i_icasigS{i}(j,:) = 0;
    % %             i_mixing_matS{i}(:,j) = 0; % according to ICA removes EEG
    % %             artifacts - just set rows of activation waveforms to be zero?
    % %
    %
    %
    %
    %             start_index = start_index+1;
    %         else
    %
    %             [f,P1] = spectralAnalysisComp(fs_data,i_icasigS{i}(j,:));
    %             [maxi,ind] = max(P1);
    %             rounded_f = round(f(ind),-1);
    %
    %             %  200 Hz frequency content, 60 Hz frequency content (added in
    %             %  120, 180 4-11-2017
    % %             if freqFilter
    % %                 if (mod(rounded_f,60) == 0 | mod(rounded_f,120) == 0 | mod(rounded_f,180) == 0)
    % %                     %if mod(rounded_f,60) == 0 || mod(rounded_f,200) == 0
    % %
    % %                     i_ica_art{i}(start_index,:) = i_icasigS{i}(j,:);
    % %                     i_ica_mix_art{i}(:,start_index) = i_mixing_matS{i}(:,j);
    % %
    % %                     i_icasigS{i}(j,:) = 0;
    % % %                     i_mixing_matS{i}(:,j) = 0;
    % %
    % %
    % %                     start_index = start_index+1;
    % %
    % %                 end
    % %             end
    %
    %         end
    %         %         if mod(rounded_f,200) == 0
    %         %             i_ica_kept{i}(start_index,:) = i_icasigS{i}(j,:);
    %         %             i_ica_mix_kept{i}(:,start_index) = i_mixing_matS{i}(:,j);
    %         %             start_index = start_index+1;
    %         %         end
    %
    %
    %     end
    
    
end

%%
%%%%%%%%%%%%%%%%%%%%%%%
% reconstruct stim artifact across channels

% make matrix of reconstruction artifacts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:numTrials
    
    recon_artifact_temp = (i_ica_mix_art{i}*i_ica_art{i})'./scale_factor;
    
    total_art(:,goods) = recon_artifact_temp;
    total_art(:,badTotal) = zeros(size(recon_artifact_temp,1),size(badTotal,2));
    
    recon_artifact_cell{i} = total_art;
    recon_artifact(:,:,i) = total_art;
    num_modes_art = size(i_ica_art{i},1);
    
end

%% subtract each one of these components

processedSig_cell = {};

processedSig = zeros(size(dataIntTime,1),size(data,2),size(numTrials,1));

total_sig = zeros(size(dataIntTime,1),size(data,2));

for i = 1:numTrials
    reconstructued_sig = ((i_mixing_matS{i}*i_icasigS{i})')./scale_factor;
    num_modes_art = size(i_ica_art{i},1);
    num_modes_kept = size(i_icasigS{i},1) - num_modes_art;
    
    % add in bad channels back
    total_sig(:,goods) = reconstructued_sig;
    total_sig(:,badTotal) = zeros(size(processedSig,1),size(badTotal,2));
    processedSig(:,:,i) = total_sig;

end

%evaluated = objective_func_ica(data,processedSig,'fs',fs_data);

evaluated = objective_func_ica_huber(data,processedSig,recon_artifact,'fs',fs_data,'plotIt',1); % try using huber loss function 

end