% ------------------------------------------------------------------------
% Get every mff files and find the awakening times from the DIN, extract EEG
% data 5 minutes before awakening. If there are DINS within 2 minutes of
% each other, use the DIN closest to awakening entry time

% Author: Thomas Vanasse
% Center for Sleep and Consciousness, University of Wisconsin - Madison
% ------------------------------------------------------------------------

eeglab;
close;

% add path to mff functions
addpath('mffs');
addpath('functions');

CHANNEL_LOCATION_FILE = 'channel_location_file/HydroCelGSN256v10.sfp';

%% input filenames from separate batch folders
inputlist = uigetfile_n_dir([pwd '/../raw_aligned_data']);

% add table name (so we aren't simultaneously i/o'ing csv file)
TABLE_name = uigetfile_n_dir(pwd,'Pick Specific DREAM REPORT FILE');

%% loop through input file list
for mff_input_file = 1:length(inputlist)
    
    TABLE = readtable(TABLE_name{1},'ReadVariableNames',true); %read experimenter input file
    filename = char(inputlist(mff_input_file)); %mff raw data filename
    
    subid = filename(end-9:end-6);
    session = filename(end-4:end-3);
    
    %% create directory if it does not already exist
    current_dir = pwd;
    subdir = [current_dir '/../sub-' subid];
    [eegdir, sesdir] = get_dirs(subdir,session);
    
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
    

    % add java file necessary for mff function
    eeglab_filepath = which('eeglab');
        
    % load absolute datetime of raw file
    aligned_file = dir(fullfile(extr_dir, '*.aligned'));
    aligned_folder = strcat(extr_dir, '/', aligned_file.name);
    
     % check if file exists, if not--mark and go to next file
    if ~isfolder(aligned_folder)
        TABLE.FILE_NOT_FOUND(de_index(1)) = 1;
        continue;
    end
    
    timing = load([aligned_folder filesep 'properties.mat']);
    absolute_datetime = timing.properties.recording_start_datetime;
    
    chaninfo = load([aligned_folder filesep 'channels.mat']);
    srate = chaninfo.channels.sampling_rate;
    
    % for each event, convert timestamp to sample
%     mff_events = mff_import_events(extr_filname);
    mff_filename = dir(fullfile(aligned_folder, '*mff_events.mat'));
    mff_events = load([aligned_folder filesep mff_filename.name]);
    mff_events = mff_events.mff_events;
    
    events_field = 'Events_DIN_1';
    
    if ~isempty(mff_events)
        events_field = 'Events_DIN_1';

        if ~isfield(mff_events, events_field)
            events_field = 'Events_Tech_Markup';
            if ~isfield(mff_events, events_field)
                events_field = 'Events_DINS';
                if ~isfield(mff_events, events_field)
                    fprintf('subject %s, %s events do not have a known field structure\n', subid, session)
                    %fieldnames(mff_events)
                end
            end
        end

        nevents = length(mff_events.(events_field).recording_timestamp);
        events = zeros(nevents, 3);  

        for ievt = 1:nevents
            datetime_start = datetime(mff_events.(events_field).recording_timestamp(ievt),...
                        'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS');
            events(ievt, 2) =  round(seconds(datetime_start - absolute_datetime) *...
                srate);
        end   

        % keep only events with code DIN100 (in theory the code for Serial Awakening)
        din100 = cellfun(@(x)(strcmp(x, 'D100')), mff_events.(events_field).code);
        din100_idx = find(din100);

         % also impose that it has at least 20 minutes and multiple DINs in one

        duration_between_dins = [events(din100_idx(1),2)/srate/60 diff(events(din100_idx,2))'/srate/60];
        timestamps = mff_events.(events_field).recording_timestamp(din100_idx)';
        timenum = zeros(1,length(timestamps));

        txt_file = [sesdir '/din_events_info.txt'];
        fileID = fopen(txt_file, 'w');
        for its = 1:length(timestamps)
            fprintf(fileID, 'dinorder-%d %s (%s %2.2f minutes)\n', ...
                its, ...
                timestamps{its}(12:16), ...
                mff_events.(events_field).code{din100_idx(its)}, ...
                duration_between_dins(its));
        end
        fclose(fileID);
        
        %%  match one din event to one awakening (closest within +/- 2 minutes)
        din_event_match = zeros(length(timestamps),4); % create array to mark
        din_event_match(:,4) = 1; % keeps track if this DIN event has the shortest time away from data entry event
        
        for time_idx = 1:length(timestamps) % for every din_event
            
                for awake_idx = 1:length(de_index) % check to see if it matches to an awakening
                % get awakening data into datetime format
                awake_date = char(TABLE.Date_Time(de_index(awake_idx)));
 
                awake_date = [awake_date(1:6) '20' awake_date(9:16)];
                
                awake = datetime(awake_date,'InputFormat','MM/dd/yyyy HH:mm');
                % get din timestamp data into datetime format
                din_timestamp = char(timestamps(time_idx));
                din_timestamp = [din_timestamp(1:10) ' ' din_timestamp(12:23)]; 
                din_timestamp = datetime(din_timestamp,'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
                
                time_interval_in_seconds = abs(etime(datevec(awake),datevec(din_timestamp)));

                if time_interval_in_seconds < 480 % +/- 8 minutes
                    din_event_match(time_idx,1) = awake_idx; %match din to specific awakening
                    din_event_match(time_idx,2) = time_interval_in_seconds; %track error (seconds)
                    break
                else 
                    continue
                    
                end
                end
        end
        
        % count the number of DIN EVENTS without awakening match
        num_zeros = sum(din_event_match(:,1)==0);
        % find repeated elements (i.e., multiple dins for one awakening)
        din_event_match(:,3) = sum(din_event_match(:,1)==din_event_match(:,1)');
        
        %% identify DIN with the shortest time from recorded awakening
        
        for k = 1:length(din_event_match(:,1))
            if ((din_event_match(k,1) > 0) && (din_event_match(k,3) > 1)) %3rd column specifices number of Awakening repeats,...
                % 4th column is valid awakening
               TABLE.DIN_EVENTS_WITHIN_TWO_MINUTES(de_index(din_event_match(k,1))) = din_event_match(k,3);
               mark_idx = k; %keeps track of index with shortest distance away from data_entry event 
               for z = 1:din_event_match(k,3)-1
                    if (din_event_match(k+1,2) < din_event_match(mark_idx,2))                 
                        din_event_match(mark_idx,4)=0;
                        din_event_match(mark_idx,3)=0; %Prevents looping again
                        mark_idx = k+1; 
                        din_event_match(mark_idx,4)=1;
                        k=k+1;
                        din_event_match(mark_idx,3)=0; %Prevents looping again
                    else 
                        din_event_match(k+1,4)=0;
                        din_event_match(k+1,3)=0; %Prevents looping again
                        k=k+1;
                    end       
                end
            else 
                continue
            end
        end
        
        fprintf('\nExtraction Done %s, %s... \n\n', subid, session)
        
    
    end 
    
    timestamps
    din_event_match
    fprintf('Found %d of %d awakenings\n', size(nonzeros(unique(din_event_match(:,1))),1), length(de_index));
    clear timestamps din_event_match

end





