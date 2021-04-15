% ------------------------------------------------------------------------
% Concatenate cleaned data to run ica again

% Author: Thomas Vanasse
% Center for Sleep and Consciousness, University of Wisconsin - Madison
% ------------------------------------------------------------------------

addpath('other/');
addpath('functions/');

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
    
    addpath(genpath(sesdir));
    
    AMICA_DIR = [local_folder{1} '/amicaout2/'];  % good to have a separate dir for each file
    if ~exist(AMICA_DIR,'dir')
        mkdir(AMICA_DIR)
    end
    
    %load index of nrem awakenings
    load([sesdir '/nrem_index'])
    
    % change to local directory to save output
    scripts_dir = pwd;
    cd(local_folder{1})
    
    % merge cleaned data and document awakening lengths after cleaning
    cleaned_lengths = [];
    
    for x = 1:length(nrem_index)
        if x == 1
           MERGEDEEG = pop_loadset([sesdir '/awakening-' num2str(nrem_index(x)) '-nrem-posICA1-cleaned.set']);
           cleaned_lengths = [cleaned_lengths MERGEDEEG.pnts];
        else 
            EEG = pop_loadset([sesdir '/awakening-' num2str(nrem_index(x)) '-nrem-posICA1-cleaned.set']);
            MERGEDEEG = pop_mergeset(MERGEDEEG,EEG); 
            cleaned_lengths = [cleaned_lengths EEG.pnts];
        end     
    end
    
    EEG = MERGEDEEG;
    clear MERGEDEEG;
    
    
    save([sesdir '/cleaned_lengths.mat'], 'cleaned_lengths');
    

    % amica15 was failing on my (TJV) machine when manually setting threads
    % to 6, 0 specifies default number of threads
    if max_threads == 0
     runamica15(EEG.data, 'num_chans', EEG.nbchan,...
    'outdir', AMICA_DIR,...
    'num_models', 1, 'num_mix_comps', 3)
        
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

    EEG = eeg_checkset(EEG, 'ica'); % update EEG.icaact
    
    EEG = pop_saveset(EEG, 'filename', ['nrem_merged_ica2.set'], 'filepath', sesdir); 
    
    % change back to scripts directory
    if ~exist([sesdir '/amicaout2'],'dir')
        mkdir([sesdir '/amicaout2'])
    end
    copyfile(AMICA_DIR, [sesdir '/amicaout2'])
    cd(scripts_dir)
end

for mff_input_file = 1:length(inputlist)

    filename = char(inputlist(mff_input_file)); %mff raw data filename
    fprintf([filename ' second ica finished\n']);
    
end

