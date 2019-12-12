% ------------------------------------------------------------------------
% Concatenate cleaned data to run ica again

% Author: Thomas Vanasse
% Center for Sleep and Consciousness, University of Wisconsin - Madison
% ------------------------------------------------------------------------

addpath('other/');

%nrem_awakenings = EEG_all.pnts/(EEG_all.srate*60*5 - EEG_all.srate*2 + 1)
%save('nrem_index.mat','nrem_index');
log_file = 'log/log_file.txt';
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
    
    % merge cleaned data
    cleaned_lengths = [];
    for x = 1:length(nrem_index)
        if x == 1
           MERGEDEEG = pop_loadset([sesdir '/awakening-' num2str(nrem_index(x)) '-cleaned_nrem.set']);
           cleaned_lengths = [cleaned_lengths MERGEDEEG.pnts];
        else 
            EEG = pop_loadset([sesdir '/awakening-' num2str(nrem_index(x)) '-cleaned_nrem.set']);
            MERGEDEEG = pop_mergeset(MERGEDEEG,EEG); 
            cleaned_lengths = [cleaned_lengths EEG.pnts];
        end     
    end
    
    EEG = MERGEDEEG;
    
    save([sesdir '/cleaned_lengths.mat'], 'cleaned_lengths');
    
    AMICA_DIR = [sesdir '/amicaout2/']; % good to have a separate dir for each file
    if ~exist(AMICA_DIR,'dir')
        mkdir(AMICA_DIR)
    end

    runamica15(EEG.data, 'num_chans', EEG.nbchan,...
        'outdir', AMICA_DIR,...
        'num_models', 1, 'num_mix_comps', 3, 'max_threads',6);

    EEG.etc.amica  = loadmodout15(AMICA_DIR);
    EEG.etc.amica.S = EEG.etc.amica.S(1:EEG.etc.amica.num_pcs, :); % Weirdly, I saw size(S,1) be larger than rank. This process does not hurt anyway.
    EEG.icaweights = EEG.etc.amica.W; % unmixing weights
    EEG.icasphere  = EEG.etc.amica.S; % sphering matrix
    EEG.icawinv = EEG.etc.amica.A; % model component matrices

    EEG = eeg_checkset(EEG, 'ica'); % update EEG.icaact
    
    EEG = pop_saveset(EEG, 'filename', ['nrem_awakening_eeg_hp_trim_merged_nobadch_ica2.set'], 'filepath', sesdir); % '_ica'
    
    
end


