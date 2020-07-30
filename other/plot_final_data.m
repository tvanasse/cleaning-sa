% ------------------------------------------------------------------------
% Get ICA cleaned data, split back to separate awakenings, then interpolate
% bad channels, and average reference data.

% Author: Thomas Vanasse
% Center for Sleep and Consciousness, University of Wisconsin - Madison
% ------------------------------------------------------------------------

addpath('other/');
addpath(genpath('functions/'));
eeglab;
close;

%% input filenames
inputlist = get_ses_dirs();

%% loop through input file list
for mff_input_file = 1:length(inputlist)
    
    sesdir = char(inputlist(mff_input_file));
   
    
    %% continue with script
    load([sesdir '/nrem_index'])
    
    EEG_all = pop_loadset([sesdir '/nrem_merged_ica2_subcomps.set']);
    
    % remove additionally identified bad channels
%     EEG = epi_log(@pop_select, EEG, 'nochannel', EEG.badchannels);
    
    %50 Hz low-pass filter
    EEG_all = pop_eegfiltnew(EEG_all, [], 50, [], 0, [], 0); 
    
    load([sesdir '/cleaned_lengths.mat']);
    cumul_sample = 1; %keep track of sample 
    
for awak = 1:length(nrem_index)
    
    %get entire five minutes
    start_sample = cumul_sample;
    end_sample = cumul_sample + cleaned_lengths(awak) - 1; 
    fprintf(['start: ' num2str(start_sample) ' end: ' num2str(end_sample) '\n'])
    cumul_sample = cumul_sample + cleaned_lengths(awak);
    
    fprintf('\nStart sample: %d\n', start_sample/500);
    fprintf('\nEnd sample: %d\n', end_sample/500);

    %% Plot topo plots before ica
    EEG = pop_select(EEG_avgref_beforeica2, 'point', [start_sample, end_sample]);

    EEG = pop_importdata('dataformat','array','data',EEG.data,...
        'srate',500,'xmin',0,'nbchan',EEG_all.nbchan, 'chanlocs', EEG_all.chanlocs);

    % load bad sections
    
    if isfile([sesdir '/awakening-' num2str(nrem_index(awak)) '-badsections.mat'])
        load([sesdir '/awakening-' num2str(nrem_index(awak)) '-badsections.mat']);
        EEG.badsections = badsections;
    else
        EEG.badsections = []; % no bad sections
    end
    

    % identify and interpoloate bad channels
    load([sesdir '/chanlocs_185.mat'])
    EEG.urchanlocs = origEEGchanlocs;   
    EEG = epi_log(@eeg_interp, EEG, EEG.urchanlocs);
    
    % average reference
    EEG = epi_log(@pop_reref, EEG, []);

    EEG = csc_eeg_plotter(EEG);
    
    waitfor(EEG)
    
end
end
