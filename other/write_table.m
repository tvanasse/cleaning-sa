% ------------------------------------------------------------------------
% Find DIN events that match to experimenter reported awakening. If there
% are multiple DIN events for one awakening, then find the closest DIN
% event to the reported one. For each awakening found, extract five minutes
% of data before awakening. Concatenate those files that were awoke in NREM
% sleep for ICA cleaning. 

% Author: Thomas Vanasse
% Center for Sleep and Consciousness, University of Wisconsin - Madison
% ------------------------------------------------------------------------

eeglab;
close;

% add functions from folder
addpath('functions');

CHANNEL_LOCATION_FILE = 'channel_location_file/HydroCelGSN256v10.sfp';

%% input filenames from separate batch folders
inputlist = uigetfile_n_dir('/Volumes/NCCAM/NCCAM/NCCAM3/SerialAwakenings/FINALS');

% add table name (so we aren't simultaneously i/o'ing csv file)
TABLE_name = uigetfile_n_dir(pwd,'Pick Specific DREAM REPORT FILE');

%% loop through input file list
for mff_input_file = 1:length(inputlist)
    
    TABLE = readtable(TABLE_name{1},'ReadVariableNames',true); %read experimenter input file
    filename = char(inputlist(mff_input_file)); %mff raw data filename
    
    subid = filename(end-9:end-6);
    session = filename(end-4:end-3);
    
    %% hard code for certain folders
%     subid = '2034'
%     session = 'T2'
    
    %% create directory if it does not already exist
    current_dir = pwd;
    subdir = [current_dir '/../sub-' subid];
    [eegdir, sesdir, timepoint] = get_dirs(subdir,session);
    
    fprintf('\nExtracting data for subject %s, %s... \n\n', subid, session)
   
    extr_dir = filename;

    %% find data entry indexes for the appropriate subject & timepoint
    de_index = [];
    for search = 1:length(TABLE.ID) 
        if ((TABLE.ID(search)==str2double(subid)) && (TABLE.Visit(search)==timepoint))
            de_index = [de_index search];
        else
        end
    end 
        
    % load absolute datetime of raw file
    aligned_file = dir(fullfile(extr_dir, '*.aligned'));
    aligned_folder = strcat(extr_dir, '/', aligned_file.name);
    
     % check if file exists, if not--mark and go to next file
    if ~isfolder(aligned_folder)
        TABLE.FILE_NOT_FOUND(de_index(1)) = 1;
        continue;
    end
    
    %% save table
    writetable(TABLE(de_index,:),[sesdir '/awakening_table.txt'])
    
end 





