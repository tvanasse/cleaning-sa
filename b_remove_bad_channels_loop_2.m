% ------------------------------------------------------------------------
% Manually remove bad channels

% Author: Thomas Vanasse
% Center for Sleep and Consciousness, University of Wisconsin - Madison
% ------------------------------------------------------------------------

CHANNEL_LOCATION_FILE = 'channel_location_file/HydroCelGSN256v10.sfp';

addpath('functions');
eeglab;
close;

%% input filenames
inputlist = get_ses_dirs();

%% loop through input file list
for mff_input_file = 1:length(inputlist)
    
    sesdir = char(inputlist(mff_input_file));
    
    EEG = pop_loadset([sesdir '/nrem_awakening_eeg_hp_trim_merged.set']);
    
    % store original channel locations (185)
    origEEGchanlocs = EEG.chanlocs;
    save([sesdir '/chanlocs_185.mat'],'origEEGchanlocs');

    EEG = csc_eeg_plotter(EEG);
    
    % wait to close eeg
    waitfor(EEG);

    EEG.badchannels = sort(EEG.csc_hidden_channels);
    
    if length(EEG.csc_hidden_channels) > 0
        figure; topoplot([],EEG.chanlocs([EEG.badchannels]), 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo); % which channels removed
        saveas(gcf, [sesdir '/nrem_badchannels'], 'tif');
        close all;      
    end
    
    % remove the bad channels, saving original channel info for later
    EEG = epi_log(@pop_select, EEG, 'nochannel', EEG.badchannels);

    EEG.urchanlocs = origEEGchanlocs;

    % save to disk.
    EEG = pop_saveset(EEG, 'filename', 'nrem_awakening_eeg_hp_trim_merged_nobadch','filepath', sesdir, 'check', 'on');

end 




