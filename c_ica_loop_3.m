% ------------------------------------------------------------------------
% Perform ICA (runica) in loop for concatenated nrem awakenings

% Author: Thomas Vanasse
% Center for Sleep and Consciousness, University of Wisconsin - Madison
% ------------------------------------------------------------------------

CHANNEL_LOCATION_FILE = 'channel_location_file/HydroCelGSN256v10.sfp';

addpath('functions');

eeglab;
close;

%% input filenames
more_files = 'Yes';
first_iter = 0;
while strcmp(more_files, 'No') ~= 1
batch_folder = uigetfile_n_dir();
if first_iter == 0;
    inputlist = uigetfile_n_dir(batch_folder);
else
    inputlist = [inputlist uigetfile_n_dir(batch_folder)];
end

more_files = questdlg('Add more subjects from different batch folder?', ...
	'??', ...
	'Yes','No','No');
first_iter = first_iter + 1;
end

% set local folder to run ica
local_folder = uigetfile_n_dir('~');

% manually set max threads - TJV: 3, TA: 6
max_threads = inputdlg('Input Max Threads',...
         'Sample', [1 50]);
max_threads = str2num(max_threads{1});


%% loop through input file list

for mff_input_file = 1:length(inputlist)
    
    TABLE = readtable('NCCAM3_06_SADreamReports_10-20-18.csv','ReadVariableNames',true); %read experimenter input file
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
    
    
    EEG = pop_loadset([sesdir '/nrem_awakening_eeg_hp_trim_merged_nobadch.set']);

    AMICA_DIR = [local_folder{1} '/amicaout/'];  % good to have a separate dir for each file
    if ~exist(AMICA_DIR,'dir')
        mkdir(AMICA_DIR)
    end
    
    % change to local directory to save output
    scripts_dir = pwd;
    cd(local_folder{1})
    
    % amica15 was failing on my (TJV) machine when manually setting threads
    % to 6
    if max_threads == 0
     runamica15(EEG.data, 'num_chans', EEG.nbchan,...
    'outdir', AMICA_DIR,...
    'num_models', 1, 'num_mix_comps', 3)
   % 'blk_min', 64, 'blk_step', 64, 'blk_max', 256); 
        
    else 
    runamica15(EEG.data, 'num_chans', EEG.nbchan,...
        'outdir', AMICA_DIR,...
        'num_models', 1, 'num_mix_comps', 3, 'max_threads', max_threads);
    end

    EEG.etc.amica  = loadmodout15(AMICA_DIR);
    EEG.etc.amica.S = EEG.etc.amica.S(1:EEG.etc.amica.num_pcs, :); 
    EEG.icaweights = EEG.etc.amica.W; % unmixing weights
    EEG.icasphere  = EEG.etc.amica.S; % sphering matrix
    EEG.icawinv = EEG.etc.amica.A; % model component matrices
    
    %EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);

    EEG = eeg_checkset(EEG, 'ica'); % update EEG.icaact
    
    EEG = pop_saveset(EEG, 'filename', [EEG.filename(1:end-4), '_ica.set'], 'filepath', sesdir); % '_ica'

    % copy files to sesdir and change back to scripts directory
    if ~exist([sesdir '/amicaout'],'dir')
        mkdir(AMICA_DIR)
    end
    copyfile(AMICA_DIR, sesdir)
    
    cd(scripts_dir)
    
    
end 




