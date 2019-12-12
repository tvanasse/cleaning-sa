%%written by Brady Riedner 07/23/2012 adapted from scripts by some one
%%better at code writing
%%% University of wisconsin
% currently assumes all bad channels are the same - just reads all segments
function bad_channels = mff_import_badchannels(meta_file)
info_file     = [meta_file filesep 'info1.xml'];
id            = fopen(info_file, 'r');
bad_channels_info = [];

while ~feof(id)
    line = fgetl(id);
    % checks if it is a section start
    if ~isempty(strfind(line, '<channels exclusion="badChannels">'))
%     di = strsplit(line,'>');
%     bad_channels_text = di{2}; can work in 2014
    [~,di] = strtok(line,'>');
    bad_channels_text = strtrim(di(2:end));
    bad_channels_info = sscanf(bad_channels_text,'%d');
    end
end
fclose(id);

categories_file     = [meta_file filesep 'categories.xml'];
id                  = fopen(categories_file, 'r');
bad_channels_cat    = [];
if exist(categories_file,'file')
    while ~feof(id)
        line = fgetl(id);
        % checks if it is a section start
        if ~isempty(strfind(line, '<channels exclusion="badChannels">'))
            [~,di] = strtok(line,'>');
            bad_channels_text = strtrim(di(2:end));
            bad_channels_cat = cat(2,bad_channels_cat,sscanf(bad_channels_text,'%d'));
        end
    end

fclose(id);
end
bad_channels = unique([bad_channels_info bad_channels_cat]);
