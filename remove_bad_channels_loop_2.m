% ------------------------------------------------------------------------
% Manually remove bad channels

% Author: Thomas Vanasse
% Center for Sleep and Consciousness, University of Wisconsin - Madison
% ------------------------------------------------------------------------

CHANNEL_LOCATION_FILE = 'channel_location_file/HydroCelGSN256v10.sfp';

%% read allfilenames

inputlist = uigetfile_n_dir;

%% loop through input file list

for mff_input_file = 1:length(inputlist)
    
    TABLE = readtable('NCCAM3_06_SADreamReports_10-20-18.csv','ReadVariableNames',true); %read experimenter input file
    filename = char(inputlist(mff_input_file)); %mff raw data filename
    subid = filename(end-9:end-6); %subject ID
    
    %% create directory if it does not already exist
    
    subdir = ['/Volumes/data/NCCAM3/SA/wDreamReport/aligned/extraction_TJV' '/sub-' subid];
    if ~exist(subdir, 'dir')
       mkdir(subdir);
       eegdir = [subdir '/eeg'];
       mkdir(eegdir);
    else 
       eegdir = [subdir '/eeg'];
    end

    session = filename(end-4:end-3);
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




