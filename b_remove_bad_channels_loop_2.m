% ------------------------------------------------------------------------
% Manually remove bad channels

% Author: Thomas Vanasse
% Center for Sleep and Consciousness, University of Wisconsin - Madison
% ------------------------------------------------------------------------

CHANNEL_LOCATION_FILE = 'channel_location_file/HydroCelGSN256v10.sfp';

addpath('functions');
eeglab;
close;

%% read allfilenames
batch_folder = uigetfile_n_dir();
inputlist = uigetfile_n_dir(batch_folder);

%% loop through input file list

for mff_input_file = 1:length(inputlist)
    
    filename = char(inputlist(mff_input_file)); %mff raw data filename
    
    k = strfind(filename, 'NCCAM_'); % get start index for filename (from full path)
    subid = filename(k+6:k+9);
    session = filename(k+11:k+12);
     
    %% create directory if it does not already exist
    current_dir = pwd;
    subdir = [current_dir '/../sub-' subid];
    
    if ~exist(subdir, 'dir')
       mkdir(subdir);
       eegdir = [subdir '/eeg'];
       mkdir(eegdir);
    else 
       eegdir = [subdir '/eeg'];
    end

    
    if strcmp(session,'T1')
        sesdir = [eegdir '/ses-1'];
        timepoint=1;
        if ~exist(sesdir)
            mkdir(sesdir);
        end
    elseif strcmp(session,'T2')
        sesdir = [eegdir '/ses-2'];
        timepoint=2;
        if ~exist(sesdir)
            mkdir(sesdir);
        end
    elseif strcmp(session,'T3')
        sesdir = [eegdir '/ses-3'];
        timepoint=3;
        if ~exist(sesdir)
            mkdir(sesdir);
        end
    else
        fprintf('session directory error\n');
        sesdir = [eegdir '/error'];
    end 
    
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




