% ------------------------------------------------------------------------
% Plot ICA activations and remove bad stretches

% Author: Thomas Vanasse
% Center for Sleep and Consciousness, University of Wisconsin - Madison
% ------------------------------------------------------------------------

addpath('other/');
addpath('functions/');

eeglab;
close;

%% input filenames
inputlist = get_ses_dirs();

%% loop through input file list
for mff_input_file = 1:length(inputlist)
    
    sesdir = char(inputlist(mff_input_file));
    
    load([sesdir '/nrem_index'])
    
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
    
    EEG_beforeclean = pop_loadset([sesdir '/nrem_awakening_eeg_hp_trim_merged_nobadch_ica.set']);
  
    if isempty(EEG_beforeclean.icaact) 
        EEG_beforeclean.icaact = (EEG_beforeclean.icaweights*EEG_beforeclean.icasphere)*EEG_beforeclean.data(EEG.icachansind,:);
    end
    
    eegplot(EEG_beforeclean.icaact,...
        'srate', EEG_beforeclean.srate,...
        'eloc_file', EEG_beforeclean.chanlocs,...
        'spacing', 5,...
        'winlength', 30,...
        'dispchans', 60,...
        'command', 'disp(''No data rejected. Use pop_select for this.'')',...
        'butlabel', 'MARK', ...
        'events', EEG_beforeclean.event);


   waitfor(gcf) %wait for epochs to be clean

   

   close all % close eeglab
        

    
end

