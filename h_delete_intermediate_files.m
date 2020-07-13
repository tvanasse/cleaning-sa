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
    
    %% remove intermediate files
    delete([sesdir '/*_eeg.set']);
    delete([sesdir '/*_eeg_hp_trim.set'])
    delete([sesdir '/nrem_awakening_eeg_hp_trim_merged.set'])
    delete([sesdir '/nrem_awakening_eeg_hp_trim_merged_nobadch.set'])
    
    delete([sesdir '/*_eeg.fdt']);
    delete([sesdir '/*_eeg_hp_trim.fdt'])
    delete([sesdir '/nrem_awakening_eeg_hp_trim_merged.fdt'])
    delete([sesdir '/nrem_awakening_eeg_hp_trim_merged_nobadch.fdt'])
    
    %% raw data
    % rmdir(filename,'s')
    
end
