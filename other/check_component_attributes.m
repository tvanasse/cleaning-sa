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
    
    pop_viewprops(EEG_avgref, 0,[1:size(EEG.icaweights,1)]);
    csc_eeg_plotter(EEG);
    
    waitfor(gcf)
    
end 




