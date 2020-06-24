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
inputlist = get_ses_dirs();

% set local folder to run ica
local_folder = uigetfile_n_dir('~');

% manually set max threads - TJV: 3, TA: 6
max_threads = inputdlg('Input Max Threads',...
         'Sample', [1 50]);
max_threads = str2num(max_threads{1});

%% loop through input file list
for mff_input_file = 1:length(inputlist)
    
    sesdir = char(inputlist(mff_input_file));
      
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
    
    
    %% if script fails to save run lines up to cd(scripts_dir)
    EEG = pop_saveset(EEG, 'filename', [EEG.filename(1:end-4), '_ica.set'], 'filepath', sesdir); % '_ica'

    % copy files to sesdir and change back to scripts directory
    if ~exist([sesdir '/amicaout'],'dir')
        mkdir([sesdir '/amicaout'])
    end
    copyfile(AMICA_DIR, [sesdir '/amicaout'])
    
    cd(scripts_dir)
    
  
    
    
end 




