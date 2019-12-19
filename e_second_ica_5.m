% ------------------------------------------------------------------------
% Concatenate cleaned data to run ica again

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
    
    %TABLE = readtable('NCCAM3_06_SADreamReports_10-20-18.csv','ReadVariableNames',true); %read experimenter input file
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
    
    EEG = pop_loadset([sesdir '/nrem_awakening_eeg_hp_trim_merged_icaprep.set']);
    
    AMICA_DIR = [sesdir '/amicaout2/']; % good to have a separate dir for each file
    if ~exist(AMICA_DIR,'dir')
        mkdir(AMICA_DIR)
    else 
        rmdir(AMICA_DIR, 's')
        mkdir(AMICA_DIR)
    end

%     runamica15(EEG.data, 'num_chans', EEG.nbchan,...
%         'outdir', AMICA_DIR,...
%         'num_models', 1, 'num_mix_comps', 3, 'max_threads',6);
    
    % amica15 was failing on my (TJV) machine when manually setting threads
    % to 6
    
    runamica15(EEG.data, 'num_chans', EEG.nbchan,...
        'outdir', AMICA_DIR,...
        'num_models', 1, 'num_mix_comps', 3);

    EEG.etc.amica  = loadmodout15(AMICA_DIR);
    EEG.etc.amica.S = EEG.etc.amica.S(1:EEG.etc.amica.num_pcs, :); % Weirdly, I saw size(S,1) be larger than rank. This process does not hurt anyway.
    EEG.icaweights = EEG.etc.amica.W; % unmixing weights
    EEG.icasphere  = EEG.etc.amica.S; % sphering matrix
    EEG.icawinv = EEG.etc.amica.A; % model component matrices

    EEG = eeg_checkset(EEG, 'ica'); % update EEG.icaact
    
    EEG = pop_saveset(EEG, 'filename', ['nrem_merged_ica2.set'], 'filepath', sesdir); % '_ica'
    
    
end


