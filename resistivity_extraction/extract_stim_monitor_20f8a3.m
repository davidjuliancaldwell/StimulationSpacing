% DJC - 1/10/2018
% script to extract the monitored stim output voltage for 20f8a3

dataStruct = struct('pair_21_20',struct('stim_current',[],'time_vec',[],'stim_data',[])...
    ,'pair_20_12',struct('stim_current',[],'time_vec',[],'stim_data',[])...
    ,'pair_22_19',struct('stim_current',[],'time_vec',[],'stim_data',[])...
    ,'pair_23_18',struct('stim_current',[],'time_vec',[],'stim_data',[])...
    ,'pair_28_4',struct('stim_current',[],'time_vec',[],'stim_data',[]));

% vector to loop through

pair_vec = {'pair_21_20','pair_20_12','pair_22_19','pair_23_18','pair_28_4'};
%%
% 21-20
for pair = pair_vec
    switch(char(pair))
        case 'pair_21_20'
            load('G:\My Drive\GRIDLabDavidShared\20f8a3\StimulationSpacing\StimSpacing_21_20_Wide.mat')
        case 'pair_20_12'
            load('G:\My Drive\GRIDLabDavidShared\20f8a3\StimulationSpacing\StimSpacing_20_12_Wide.mat')
        case 'pair_22_19'
            load('G:\My Drive\GRIDLabDavidShared\20f8a3\StimulationSpacing\StimSpacing_22_19_Wide.mat')
        case 'pair_23_18'
            load('G:\My Drive\GRIDLabDavidShared\20f8a3\StimulationSpacing\StimSpacing_23_18_Wide.mat')
        case 'pair_28_4'
            load('G:\My Drive\GRIDLabDavidShared\20f8a3\StimulationSpacing\StimSpacing_28_4_Wide.mat')
    end
    % get sampling rates
    fs_data = Wave.info.SamplingRateHz;
    fs_stim = Stim.info.SamplingRateHz;
    
    % stim data
    stim = Stim.data;
    
    % current data
    sing = Sing.data;
    
    %% plot stim channels if interested
    
    
    %% Sing looks like the wave to be delivered, with amplitude in uA
    
    
    % build a burst table with the timing of stimuli
    bursts = [];
    
    % first channel of current
    Sing1 = sing(:,1);
    fs_sing = Sing.info.SamplingRateHz;
    
    samplesOfPulse = round(2*fs_stim/1e3);
    
    Sing1Mask = Sing1~=0;
    dmode = diff([0 Sing1Mask' 0 ]);
    
    dmode(end-1) = dmode(end);
    
    bursts(2,:) = find(dmode==1);
    bursts(3,:) = find(dmode==-1);
    
    singEpoched = squeeze(getEpochSignal(Sing1,(bursts(2,:)-1),(bursts(3,:))+1));
    t = (0:size(singEpoched,1)-1)/fs_sing;
    t = t*1e3;
    
    %% Plot stims with info from above, and find the delay!
    
    stim1stChan = stim(:,1);

    % get the delay in stim times - looks to be 7 samples or so
    delay = round(0.2867*fs_stim/1e3);
   
    % plot the appropriately delayed signal
    stimTimesBegin = bursts(2,:)-1+delay - 20;
    stimTimesEnd = bursts(3,:)-1+delay+120;
    stim1Epoched = squeeze(getEpochSignal(stim1stChan,stimTimesBegin,stimTimesEnd));
    
    t = (-20:size(stim1Epoched,1)-20-1)/fs_stim;
    t = t*1e3;
    stimdata = stim1Epoched;
    time_vec = t;
    stim_current = max(max(Sing.data));
    
    switch(char(pair))
        case 'pair_21_20'
            dataStruct.pair_21_20.stim_current = stim_current;
            dataStruct.pair_21_20.time_vec = time_vec;
            dataStruct.pair_21_20.stim_data = stimdata;
        case 'pair_20_12'
            dataStruct.pair_20_12.stim_current = stim_current;
            dataStruct.pair_20_12.time_vec = time_vec;
            dataStruct.pair_20_12.stim_data = stimdata;
        case 'pair_22_19'
            dataStruct.pair_22_19.stim_current = stim_current;
            dataStruct.pair_22_19.time_vec = time_vec;
            dataStruct.pair_22_19.stim_data = stimdata;
        case 'pair_23_18'
            dataStruct.pair_23_18.stim_current = stim_current;
            dataStruct.pair_23_18.time_vec = time_vec;
            dataStruct.pair_23_18.stim_data = stimdata;
        case 'pair_28_4'
            dataStruct.pair_28_4.stim_current = stim_current;
            dataStruct.pair_28_4.time_vec = time_vec;
            dataStruct.pair_28_4.stim_data = stimdata;
    end
    
    
end

save('20f8a3_stimOutput','dataStruct');