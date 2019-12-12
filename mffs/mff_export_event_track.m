% some modifications have been made to improve performance 12/13/14
function mff_export_event_track(meta_file, event_track_name, events)
    absolute_timestamp = mff_find_absolute_timestamp(meta_file);
    
    event_track_file    = ['Events_' strrep(event_track_name, ' ', '_') '.xml'];
    event_track_file_id = fopen([meta_file filesep event_track_file], 'W');
    
    fprintf(event_track_file_id, '<?xml version="1.0" encoding="iso-8859-1" standalone="yes"?>\n');
    fprintf(event_track_file_id, '<eventTrack xmlns="http://www.egi.com/event_mff" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">\n');
    fprintf(event_track_file_id, '\t<name>%s</name>\n', event_track_name);
    fprintf(event_track_file_id, '\t<trackType>EVNT</trackType>\n');
    
    num_events = length(events.code);
    
    for event_num = 1:num_events
        % section needs to be modified to allow easy entry of events
        % depending on desired format
%         onset    = round(1000 * events.onset(event_num));  % need to be
%         duration = round(1000000 * events.duration(event_num));
        
        fprintf(event_track_file_id, '\t<event>\n');
        %fprintf(event_track_file_id, '\t\t<beginTime>%s%s%s</beginTime>\n', datestr(addtodate(absolute_timestamp.datenum, onset, 'millisecond'), 'yyyy-mm-ddTHH:MM:SS.FFF'), absolute_timestamp.offset, absolute_timestamp.timezone);
        fprintf(event_track_file_id, '\t\t<beginTime>%s%s%s</beginTime>\n',events.recording_timestamp{event_num},events.recording_timestamp_submilliseconds{event_num},events.recording_timestamp_timezone{event_num});
        %fprintf(event_track_file_id, '\t\t<duration>%lu000</duration>\n', duration);
        fprintf(event_track_file_id, '\t\t<duration>%lu</duration>\n', events.duration(event_num));
        fprintf(event_track_file_id, '\t\t<code>%s</code>\n', events.code{event_num});
        fprintf(event_track_file_id, '\t\t<label>%s</label>\n', events.label{event_num});
        fprintf(event_track_file_id, '\t\t<description>%s</description>\n', events.description{event_num});
        fprintf(event_track_file_id, '\t\t<sourceDevice>%s</sourceDevice>\n', events.sourceDevice{event_num});
        
        num_keys = length(events.keys(event_num).keyCode);
        
        if num_keys > 0
            fprintf(event_track_file_id, '\t\t<keys>\n');
            
            for key_num = 1:num_keys
                fprintf(event_track_file_id, '\t\t\t<key>\n');
                fprintf(event_track_file_id, '\t\t\t\t<keyCode>%s</keyCode>\n', events.keys(event_num).keyCode{key_num});
                switch events.keys(event_num).dataType{key_num}
                    case 'string'
                fprintf(event_track_file_id, '\t\t\t\t<data dataType=\"%s\">%s</data>\n', events.keys(event_num).dataType{key_num}, events.keys(event_num).data{key_num});
                    case 'long'
                fprintf(event_track_file_id, '\t\t\t\t<data dataType=\"%s\">%lu</data>\n', events.keys(event_num).dataType{key_num}, events.keys(event_num).data{key_num});
                end       
                fprintf(event_track_file_id, '\t\t\t</key>\n');
            end
            
            fprintf(event_track_file_id, '\t\t</keys>\n');
        end
        
        fprintf(event_track_file_id, '\t</event>\n');
    end
    
    fprintf(event_track_file_id, '</eventTrack>\n');
    
    fclose(event_track_file_id);fclose all;
end

