% ------------------------------------------------------------------------
% Plot ICA activations and remove bad stretches

% Author: Thomas Vanasse
% Center for Sleep and Consciousness, University of Wisconsin - Madison
% ------------------------------------------------------------------------

addpath('other/');
addpath('functions/');

eeglab;
close;

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
    
    load([sesdir '/nrem_index'])
    
    EEG = pop_loadset([sesdir '/nrem_awakening_eeg_hp_trim_merged_nobadch_ica.set']);
    
    EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:); %from eeg_checkset
    
    eegplot(EEG.icaact,...
    'srate', EEG.srate,...
    'eloc_file', EEG.chanlocs,...
    'spacing', 5,...
    'winlength', 30,...
    'dispchans', 128,...
    'command', 'disp(''No data rejected. Use pop_select for this.'')',...
    'butlabel', 'MARK', ...
    'events', EEG.event);

    waitfor(gcf) %wait for epochs to be clean
    
    eegh %history function
    
    close all % close eeglab
    
    % data storing the removed sections
    EEG.badsections = TMPREJ(:,1:2);
    
    % calculate percentage of removed sections
    sum_rmv = 0;
    for rmv = 1:length(TMPREJ(:,1:2))
        sum_rmv = sum_rmv + (TMPREJ(rmv,2)-TMPREJ(rmv,1));
    end
    
    fid = fopen([sesdir '/log.txt'], 'at');
    fprintf(fid, 'Percentage of data removed: %.2f\n', sum_rmv/EEG.pnts*100);
    fclose(fid);
    
    % remove bad streches, but keep badsections data in 'badsections' field
    EEG = epi_log(@pop_select, EEG, 'nopoint', TMPREJ(:,1:2));
    
    EEG = pop_saveset(EEG, 'filename', 'nrem_awakening_eeg_hp_trim_merged_icaprep.set','filepath',sesdir);

end

