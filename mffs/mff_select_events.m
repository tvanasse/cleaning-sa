% function to adjust events based on sample start / end of file
%  recording_datetime_start = datetime(events.Events_MMN.recording_timestamp(1),'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS') - minutes(3);
%  recording_datetime_end = datetime(events.Events_MMN.recording_timestamp(2),'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS') +  seconds(events.Events_MMN.duration(2)/1000000000) + minutes(1)
%
%starttime   = datetime(recording_timestamp_start,'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS');
%endtime     = datetime(recording_timestamp_end,'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS');
function newevents = mff_select_events_and_adjust_timing(meta_file,recording_datetime_start, recording_datetime_end)

events      = mff_import_events(meta_file);
eventtracks = fieldnames(events);

for t = 1:length(eventtracks)
    
    
    eventfields  = fieldnames(events.(eventtracks{t}));
    other_event_fields  = eventfields(3:end);
    
    newevents.(eventtracks{t}).TrackName = events.(eventtracks{t}).TrackName;
    newevents.(eventtracks{t}).trackType = events.(eventtracks{t}).trackType;
    
    % copies the rest of the event info from the actual events other than
    % the TrackName and the trackType
    if ~isempty(other_event_fields) % put event times into easy to use datetime formant
        eventstarts = datetime(events.(eventtracks{t}).recording_timestamp,'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS');
        
        eventends   = eventstarts + seconds(events.(eventtracks{t}).duration/1e9);  % add durations translated into seconds
        
        % must cut the evens
        keepindex = (eventstarts > recording_datetime_start &  eventends   < recording_datetime_end);
        % add a new fields with seconds from onset
        newevents.(eventtracks{t}).onsetsecs = seconds(eventstarts(keepindex)-recording_datetime_start);
        for oefi = 1:length(other_event_fields)
            
            % if statement deals with the events that have more than one
            % element (i.e. recording_datevec)
            
            % if statement deals with the events that have more than one
            % element (i.e. recording_datevec)
            if min(size(events.(eventtracks{t}).(other_event_fields{oefi}))) == 1;
                newevents.(eventtracks{t}).(other_event_fields{oefi}) = ...
                    events.(eventtracks{t}).(other_event_fields{oefi})(keepindex); % only those two events
            else
                newevents.(eventtracks{t}).(other_event_fields{oefi}) = ...
                    events.(eventtracks{t}).(other_event_fields{oefi})(keepindex,:); % add only those two events
            end
        end
    end
end