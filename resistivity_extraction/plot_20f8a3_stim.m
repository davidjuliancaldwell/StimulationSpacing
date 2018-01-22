% DJC - 1/10/2018
% script to plot the monitored stimulation voltage from the 20f8a3 1.2 ms
% pulse experiments
%%
pair_vec = {'pair_21_20','pair_20_12','pair_22_19','pair_23_18','pair_28_4'};
% 21 20 pair
for pair = pair_vec
    pair = char(pair);
    figure
    plot(eval(strcat('dataStruct.',pair,'.time_vec')),eval(strcat('dataStruct.',pair,'.stim_data')));
    xlabel('time (ms)');
    ylabel('Voltage (V)');
    pair_title = strrep(pair,'_','\_');
    title([pair_title ' Monitored Output Voltage for Current = ' num2str(eval(strcat('dataStruct.',pair,'.stim_current'))) ' uA'])
end

