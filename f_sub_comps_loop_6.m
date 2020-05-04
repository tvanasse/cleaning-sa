% ------------------------------------------------------------------------
% Remove bad independent components with help of ICLabel

% Author: Thomas Vanasse
% Center for Sleep and Consciousness, University of Wisconsin - Madison
% ------------------------------------------------------------------------

CHANNEL_LOCATION_FILE = 'channel_location_file/HydroCelGSN256v10.sfp';

addpath('functions');
eeglab;
close;

log_file = 'log/log_file.txt';

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
    
    
    EEG = pop_loadset([sesdir '/nrem_merged_ica2.set']);
    

    % average reference for ICLabel
    EEG_avgref = epi_log(@pop_reref, EEG, []);
    
    EEG_avgref = iclabel(EEG_avgref); % prior to computing features, each dataset was converted to a common average reference
    
    class_probs = EEG_avgref.etc.ic_classification.ICLabel.classifications;

    pct_thr = 0.35; % classification threshold
    brain_cutoff = .30; % if brain contains anything over 30% in prob, keep component
    ic_classifications = zeros(7,1);
    ic_artifacts_all = [];
    ic_good_all = [];
    for i =1:length(class_probs(:,1))
        if class_probs(i,1) > pct_thr
            ic_classifications(1,1) = ic_classifications(1,1) + 1; 
            ic_good_all = [ic_good_all i]; %keep track of index       
        elseif class_probs(i,2) > pct_thr & class_probs(i,1) < brain_cutoff
            ic_classifications(2,1) = ic_classifications(2,1) + 1; 
            ic_artifacts_all = [ic_artifacts_all i];
        elseif class_probs(i,3) > pct_thr & class_probs(i,1) < brain_cutoff
            ic_classifications(3,1) = ic_classifications(3,1) + 1; 
            ic_artifacts_all = [ic_artifacts_all i];
        elseif class_probs(i,4) > pct_thr & class_probs(i,1) < brain_cutoff
            ic_classifications(4,1) = ic_classifications(4,1) + 1; 
            ic_artifacts_all = [ic_artifacts_all i];
        elseif class_probs(i,5) > pct_thr & class_probs(i,1) < brain_cutoff
            ic_classifications(5,1) = ic_classifications(5,1) + 1; 
            ic_artifacts_all = [ic_artifacts_all i];
        elseif class_probs(i,6) > pct_thr & class_probs(i,1) < brain_cutoff
            ic_classifications(6,1) = ic_classifications(6,1) + 1; 
            ic_artifacts_all = [ic_artifacts_all i];
        elseif class_probs(i,7) > pct_thr & class_probs(i,1) < brain_cutoff
            ic_classifications(7,1) = ic_classifications(7,1) + 1; 
            ic_artifacts_all = [ic_artifacts_all i];
        else 
            ic_good_all = [ic_good_all i]; %keep track of index   
        end 
    end
    
    % plot bar plot
    X = categorical(EEG_avgref.etc.ic_classification.ICLabel.classes);
    bar(X,ic_classifications)
    saveas(gcf,fullfile(sesdir,'ic_labels_2'), 'png');
    close(gcf)
    
    pop_viewprops(EEG_avgref, 0,ic_good_all);
    pop_eegplot(EEG,0,0,0);
    
    waitfor(gcf)
    
    manual_remove_components = [];
    manual_remove_components = inputdlg('REMOVE THESE COMPONENTS',...
             'Sample', [1 50]);
    manual_remove_components = str2num(manual_remove_components{1});
    
    if ~isempty(manual_remove_components)
        ic_artifacts_all = [ic_artifacts_all manual_remove_components]; %add to artifiacts
    end
    
    pop_viewprops(EEG_avgref, 0 ,ic_artifacts_all);
    pop_eegplot(EEG,0,0,0);
    waitfor(gcf)
    
    manual_keep_components = [];
    manual_keep_components = inputdlg('KEEP THESE COMPONENTS',...
             'Sample', [1 50]);
    manual_keep_components = str2num(manual_keep_components{1});
    
    if ~isempty(manual_keep_components)
        for x = 1:length(manual_keep_components)
            ids = find(ic_artifacts_all==manual_keep_components(x));
            ic_artifacts_all(ids) = [];
        end
    end
    
    % subtract artifactual components (from non-ref-averaged EEG)
    save([sesdir '/ic_artifacts.mat'],'ic_artifacts_all');
    EEG_subcomps = pop_subcomp(EEG, ic_artifacts_all, 1);
    pop_saveset(EEG_subcomps, 'filename',[EEG.filename(1:end-4) '_subcomps.set'],'filepath', sesdir);

    
end 




