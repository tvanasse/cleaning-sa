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
inputlist = get_ses_dirs();

%% loop through input file list
for mff_input_file = 1:length(inputlist)
    
    sesdir = char(inputlist(mff_input_file));
    
    EEG = pop_loadset([sesdir '/nrem_merged_ica2.set']);
    
    if isempty(EEG.icaact) 
        EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
    end
    
    eegplot(EEG.icaact,...
    'srate', EEG.srate,...
    'eloc_file', EEG.chanlocs,...
    'spacing', 5,...
    'winlength', 30,...
    'dispchans', 60,...
    'command', 'disp(''No data rejected. Use pop_select for this.'')',...
    'butlabel', 'MARK', ...
    'events', EEG.event);
    
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
    %saveas(gcf,fullfile(sesdir,'ic_labels_2'), 'png');
    close(gcf)
    
    pop_viewprops(EEG_avgref, 0,ic_good_all);
    csc_eeg_plotter(EEG);
    
    waitfor(gcf)
    
    manual_remove_components = [];
    manual_remove_components = inputdlg('REMOVE THESE COMPONENTS',...
             'Sample', [1 50]);
    manual_remove_components = str2num(manual_remove_components{1});
    
    if ~isempty(manual_remove_components)
        ic_artifacts_all = [ic_artifacts_all manual_remove_components]; %add to artifiacts
    end
    
    pop_viewprops(EEG_avgref, 0 ,ic_artifacts_all);
    
    eegplot(EEG.icaact,...
    'srate', EEG.srate,...
    'eloc_file', EEG.chanlocs,...
    'spacing', 5,...
    'winlength', 30,...
    'dispchans', 60,...
    'command', 'disp(''No data rejected. Use pop_select for this.'')',...
    'butlabel', 'MARK', ...
    'events', EEG.event);

    csc_eeg_plotter(EEG);
    
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
    % save([sesdir '/ic_artifacts.mat'],'ic_artifacts_all');
    % EEG_subcomps = pop_subcomp(EEG, ic_artifacts_all, 1);
    % pop_saveset(EEG_subcomps, 'filename',[EEG.filename(1:end-4) '_subcomps.set'],'filepath', sesdir);

    
end 




