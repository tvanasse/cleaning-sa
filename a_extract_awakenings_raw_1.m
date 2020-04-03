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

%% input filenames
inputlist = uigetfile_n_dir;

%% loop through input file list
for mff_input_file = 1:length(inputlist)
    
    TABLE = readtable('NCCAM3_06_SADreamReports_10-20-18.csv','ReadVariableNames',true); %read experimenter input file
    filename = char(inputlist(mff_input_file)); %mff raw data filename
    
    k = strfind(filename, 'NCCAM_'); % get start index for filename (from full path)
    subid = filename(k+6:k+9);
    session = filename(k+11:k+12);
    
    %% create directory if it does not already exist
    subdir = ['/Volumes/data-2/NCCAM3/SA/wDreamReport/aligned/extraction_TJV' '/sub-' subid];
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
        else 
            f = warndlg('Subject Folder Already Exists','Warning');
        end
    elseif strcmp(session,'T2')
        sesdir = [eegdir '/ses-2'];
        timepoint=2;
        if ~exist(sesdir)
            mkdir(sesdir);
        else 
            f = warndlg('Subject Folder Already Exists','Warning');
        end
    elseif strcmp(session,'T3')
        sesdir = [eegdir '/ses-3'];
        timepoint=3;
        if ~exist(sesdir)
            mkdir(sesdir);
        else 
            f = warndlg('Subject Folder Already Exists','Warning');
        end
    else
        fprintf('session directory error\n');
        sesdir = [eegdir '/error'];
    end 
    
    fprintf('\nExtracting data for subject %s, %s... \n\n', subid, session)

    mff_name = dir(fullfile(filename, '*.mff'));
    extr_filname = strcat(filename, '/', mff_name.name);
    k = strfind(filename, 'NCCAM_'); % get start index for filename (from full path)
   
    %% manually set mff file (if data is split into two parts)
%     mff_file = uigetfile_n_dir;
%     extr_filname = char(mff_file(1));
   
    extr_dir = filename;

    %% find data entry indexes for the appropriate subject & timepoint
    de_index = [];
    for search = 1:length(TABLE.ID) 
        if ((TABLE.ID(search)==str2double(subid)) && (TABLE.Visit(search)==timepoint))
            de_index = [de_index search];
        else
        end
    end 
    
    % check if file exists, if not--mark and go to next file
    if ~isfolder(extr_filname)
        TABLE.FILE_NOT_FOUND(de_index(1)) = 1;
        continue;
    end

    % add java file necessary for mff function
    eeglab_filepath = which('eeglab');
    javaaddpath([eeglab_filepath(1:end-8) 'plugins/MFFimport2.3/mffimport/MFF-1.2.jar']);
            
    % load absolute datetime of raw file
    aligned_file = dir(fullfile(extr_dir, '*.aligned'));
    aligned_folder = strcat(extr_dir, '/', aligned_file.name);
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
                    
                [Y, M, filename, H, MN, S] = datevec(awake - din_timestamp);
                seconds_diff = abs(H*3600 + MN*60 + S); 

                if seconds_diff < 300 % +/- 5 minutes
                    din_event_match(time_idx,1) = awake_idx; %match din to specific awakening
                    din_event_match(time_idx,2) = seconds_diff; %track error (seconds)
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
        
        % set din_events (without multiple per awakening) to one
        for k = 1:length(din_event_match(:,1))
            if ((din_event_match(k,1) > 0) && (din_event_match(k,3) == 1 ))
                TABLE.DIN_EVENTS_WITHIN_TWO_MINUTES(de_index(din_event_match(k,1))) = 1; 
            end 
        end 
        
        %% save five minutes of data for each awakening
        
        ent_matched_awakening = 1; %create entry matched awakening number for file names
        n2n3_iter = 0; %for merging nrem
        nrem_index = [];
        for i = 1:length(timestamps)
            %if events(din100_idx(i),2)-(srate*60*5) > 0 %check if event did not occur 5 min. to start time, skip otherwise
            if ((din_event_match(i,1) > 0) && din_event_match(i,4) > 0) %DIN matches with awakening AND is the shortest time away from recoreded awakening   
                fprintf("DIN %d \n");
                %mark time away from 
                TABLE.DIFF_FROM_EVENT_SECONDS(de_index(din_event_match(i,1))) = din_event_match(i,2);
                TABLE.ENTRY_MATCHED_AWAKENING_NO(de_index(din_event_match(i,1))) = ent_matched_awakening;
                     
                %EXTRACTING 5 MINUTES BEFORE AWAKENING
                samples_before_awakening = mff_read_samples(extr_filname,'all',...
                    events(din100_idx(i),2)-(srate*60*5), ...
                    events(din100_idx(i),2)-1);

                aligned_file = dir(fullfile(extr_dir, '*.aligned'));
                aligned_folder = strcat(extr_dir, '/', aligned_file.name);
                aligned_events_raw = load([aligned_folder filesep 'alignedevents.mat']);
                egi_offset = aligned_events_raw.alignedevents.egi_start_offset;

%                 the following code would read .raw files instead,
%                 accounting for egi offset              
%                 % get mffraw files
                f = dir(fullfile(aligned_folder,'MFF*'));
                C = struct2cell(f);
                input = {};
                for chan = 1:257
                    input{end+1} = C{1,chan};
                end 
                samples_before_awakneing = readalignmentraw(aligned_folder, input, events(din100_idx(i),2)-(srate*60*5), events(din100_idx(i),2)-1);               
                
                eloc = readlocs(CHANNEL_LOCATION_FILE);

                %eloc_no257 = eloc(1:257); %remove channel 257... so it won't affect bad channel removal
                EEG = pop_importdata('dataformat','array','data',samples_before_awakening,...
                    'srate',500,'xmin',0,'nbchan',257,'setname', sprintf('awakening-%d',ent_matched_awakening), 'chanlocs', eloc);

                EEG.condition = sprintf('session_start_time: %s, awakening_timestamp: %s, DIN_code: %s, time_since_last_event (min): %2.2f)', ...
                    absolute_datetime, ... %start time of scan session (night)
                    timestamps{i}(12:16), ... 
                    mff_events.(events_field).code{din100_idx(i)}, ...
                    duration_between_dins(i));
                EEG.subject = subid;
                
                % convert 256 to "inseide 185 channels
                load('channel_location_file/inside185ch.mat');
                EEG.chanlocs = EEG.chanlocs(inside185ch);
                EEG.nbchan = 185;
                EEG.data = EEG.data(inside185ch,:);

                EEG = pop_saveset(EEG,sprintf('awakening-%d_eeg',ent_matched_awakening), sesdir);                 
                
                if ~isfile([aligned_folder '/alignedscoring.raw'])
                    TABLE.N1(de_index(din_event_match(i,1))) = 'NO ALIGNED SCORING FILE FOUND';
                    
                else 
                    scoring = readalignmentraw(aligned_folder, {'alignedscoring.raw'}, ...
                        events(din100_idx(i),2)-(srate*60*5),...
                        events(din100_idx(i),2) - 1);
                    plot(scoring)
                    saveas(gcf, [sesdir '/awakening-' num2str(ent_matched_awakening) '-scoring.png'], 'png');
                    close all;
                    
                    TABLE.N1(de_index(din_event_match(i,1))) = any(scoring(:) == -1);
                    TABLE.N2(de_index(din_event_match(i,1))) = any(scoring(:) == -2);
                    TABLE.N3(de_index(din_event_match(i,1))) = any(scoring(:) == -3);
                    TABLE.REM(de_index(din_event_match(i,1))) = any(scoring(:) == 1);
                    TABLE.WAKE(de_index(din_event_match(i,1))) = any(scoring(:) == 0);
                end
                
                writetable(TABLE,'NCCAM3_06_SADreamReports_10-20-18.csv');
                
                % high-pass Filter (1 Hz)
                EEG = pop_eegfiltnew(EEG, 1, [], [], 0, [], 0);
                 
                % clean Line Noise
                EEG = pop_cleanline(EEG,'ChanCompIndices',[1:EEG.nbchan],'SignalType','Channels','LineFrequencies',60);
                
                % trim 1 second from beginning and end
                EEG = pop_select(EEG,'point', [0 + EEG.srate*1, EEG.pnts - EEG.srate*1]); 
                pop_saveset(EEG,sprintf('awakening-%d_eeg_hp_trim',ent_matched_awakening), sesdir);
                
                % assess if awakening occured in N2/N3 (look at 30 seconds
                % before din event
                awakening_n2n3_samples = find(scoring(EEG.pnts-30*EEG.srate:EEG.pnts) == -2 | scoring(EEG.pnts-30*EEG.srate:EEG.pnts) == -3);
                    
                % five minutes cannot contain any rem
                if (~isempty(awakening_n2n3_samples) & TABLE.REM(de_index(din_event_match(i,1))) ~= 1)
                    nrem_index = [nrem_index ent_matched_awakening];
                    
                    % find N2/N3 samples only, [ignore last 30 seconds of
                    % data] because scoring may be misleading here because
                    % data is epoched every 30 seconds
                    n2n3_samples = find(scoring(1:end-30*EEG.srate) ~= -2 & scoring(1:end-30*EEG.srate) ~= -3);
                    
                    if n2n3_iter == 0
                        MERGEDEEG = EEG;
                    else
                        MERGEDEEG = pop_mergeset(MERGEDEEG,EEG);
                    end
                    n2n3_iter = n2n3_iter + 1;                 
                end   
               ent_matched_awakening = ent_matched_awakening + 1;                         
            else 
                fprintf('No Awakening Match\n');
            end          
        end
    
    save([sesdir '/nrem_index.mat'],'nrem_index');
    pop_saveset(MERGEDEEG,'filename', 'nrem_awakening_eeg_hp_trim_merged','filepath', sesdir); 
        
    else 
        
        fprintf('subject %s, %s mff file has no events to import\n', subid, session);
        TABLE.MFF_FILE_HAS_NO_EVENTS_TO_IMPORT(de_index(1)) = 1;
        writetable(TABLE,'NCCAM3_06_SADreamReports_10-20-18.csv');

    end 
    
end 





