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
    
    k = strfind(filename, 'NCCAM_'); % get start index for filename (from full path)
    subid = filename(k+6:k+9);
    session = filename(k+11:k+12);
    
    %% create directory if it does not already exist
    subdir = ['/Volumes/data/NCCAM3/SA/wDreamReport/aligned/extraction_TJV' '/sub-' subid];
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
    
    load([sesdir '/nrem_index'])
    
    EEG_all = pop_loadset([sesdir '/nrem_awakening_eeg_hp_trim_merged_nobadch_ica.set']);
    
    if isempty(EEG_all.icaact) 
        EEG_all.icaact = (EEG_all.icaweights*EEG_all.icasphere)*EEG_all.data(EEG_all.icachansind,:);
    end
    
    for awak = 1:length(nrem_index)

        % get entire five minutes
        start_sample = (awak-1)*(EEG_all.srate*60*5 - EEG_all.srate*2 + 1);
        end_sample = awak*(EEG_all.srate*60*5 - EEG_all.srate*2 + 1); %each extraction is 4 min. 58seconds

        EEG = pop_select(EEG_all, 'point', [start_sample, end_sample]);
        
        EEG = pop_importdata('dataformat','array','data',EEG.data,...
        'srate',500,'xmin',0,'nbchan',EEG_all.nbchan, 'chanlocs', EEG_all.chanlocs);
    
        % get ica weights
        EEG.icaweights = EEG_all.icaweights;
        if awak == 1
            EEG.icaact = EEG_all.icaact(:, 1:end_sample);
        else 
            EEG.icaact = EEG_all.icaact(:, start_sample:end_sample);
        end
    
        EEG.icasphere = EEG_all.icasphere;

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
        
        if size(TMPREJ,1) > 0
            % data storing the removed sections/bad channels
            EEG.badsections = TMPREJ(:,1:2);
            badsections = EEG.badsections;
            save([sesdir '/awakening-' num2str(nrem_index(awak)) '-badsections.mat'], 'badsections');

            % calculate percentage of removed sections
            sum_rmv = 0;
            for rmv = 1:size(TMPREJ,1)
                sum_rmv = sum_rmv + (TMPREJ(rmv,2)-TMPREJ(rmv,1));
                fprintf(num2str(rmv));
            end

            fid = fopen([sesdir '/log.txt'], 'at');
            fprintf(fid, 'Percentage of data removed, awakening-%d: %.2f\n', nrem_index(awak), sum_rmv/EEG.pnts*100);
            fclose(fid);
        
        % remove bad streches, but keep badsections data in 'badsections' field
        EEG = epi_log(@pop_select, EEG, 'nopoint', TMPREJ(:,1:2));
        
        end 

        EEG = pop_saveset(EEG, 'filename', sprintf('awakening-%d-nrem-posICA1-cleaned',nrem_index(awak)),'filepath',sesdir);
    
    end
    
end

