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
% TABLE_name = uigetfile_n_dir(pwd,'Pick Specific DREAM REPORT FILE');

%% loop through input file list
for mff_input_file = 1:length(inputlist)
    
%     TABLE = readtable(TABLE_name{1},'ReadVariableNames',true); %read experimenter input file
    filename = char(inputlist(mff_input_file)); %mff raw data filename
    
    subid = filename(end-9:end-6);
    session = filename(end-4:end-3);
    
    %% hard code for certain folders
%     subid = '2034'
%     session = 'T2'
    
    %% create directory if it does not already exist
    current_dir = pwd;
    subdir = [current_dir '/../../sub-' subid];
    [eegdir, sesdir, timepoint] = get_dirs(subdir,session);
    
    fprintf('\nExtracting data for subject %s, %s... \n\n', subid, session)
   
    extr_dir = filename;

    %% find data entry indexes for the appropriate subject & timepoint
    de_index = [];
%     for search = 1:length(TABLE.ID) 
%         if ((TABLE.ID(search)==str2double(subid)) && (TABLE.Visit(search)==timepoint))
%             de_index = [de_index search];
%         else
%         end
%     end 
        
    % load absolute datetime of raw file
    aligned_file = dir(fullfile(extr_dir, '*.aligned'));
    aligned_folder = strcat(extr_dir, '/', aligned_file.name);
    
  
    scoring = readalignmentraw(aligned_folder, {'alignedscoring.raw'});
    
    fprintf('%d Hours\n', length(scoring)/(500*60*60))
    fprintf('%d Samples\n', length(scoring))

    plot(scoring,'LineWidth',2);
    hold on;
%     for i = 1:length(timestamps)
%            if ((din_event_match(i,1) > 0) && din_event_match(i,4) > 0) %DIN matches with awakening AND is the shortest time away from recoreded awakening   
%                 xline(events(din100_idx(i),2),'LineWidth',3,'Color','r');
%            end
%     end
 
    
    set(gcf, 'Units', 'Inches', 'Position', [0, 0, 18, 7]);

%     saveas(gcf, [sesdir '/all_awakenings.png'], 'png');
%     close all;
    uiwait
    
    
    %% clear timestamps/din_event_match
    clear timestamps
    clear din_event_match
   
end







