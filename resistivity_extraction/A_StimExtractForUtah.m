%% DJC - 3-24-2017 - to process .mat file for UTAH and Northeastern folks 

% clear workspace
close all; clear all; clc

% set input output working directories
%CODE_DIR = fullfile(myGetenv('gridlab_dir'));
%scripts_path =  strcat(CODE_DIR,'\Experiment\BetaTriggeredStim\scripts');

% add path for scripts to work with data tanks

% subject directory, change as needed
SUB_DIR = fullfile(myGetenv('subject_dir'));


structureData = uiimport('-file');

stim_chans = input('what are the stim chans? e.g. [12 22] \n');

Sing = structureData.Sing;
Stim = structureData.Stim;
Stm0 = structureData.Stm0;
Eco1 = structureData.ECO1;
Eco2 = structureData.ECO2;
Eco3 = structureData.ECO3;
Eco4 = structureData.ECO4;

Wave.info = structureData.ECO1.info; 
Wave.data = [Eco1.data Eco2.data Eco3.data Eco4.data];

OUTPUT_DIR = 'C:\Users\djcald.CSENETID\Data\ConvertedTDTfiles\2a8d14';
save(fullfile(OUTPUT_DIR, ['stim_widePulse_',num2str(stim_chans(1)),'_',num2str(stim_chans(2))]),'Sing','Stim','Stm0','Wave');

