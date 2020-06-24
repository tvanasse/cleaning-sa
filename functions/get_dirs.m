function [eegdir,sesdir] = get_dirs(subdir,session)
%get_dirs gets session specific and subject specific directorys, if they
%don't exist, then they are created
%   Detailed explanation goes here
 
if ~exist(subdir, 'dir')
   mkdir(subdir);
   eegdir = [subdir '/eeg'];
   mkdir(eegdir);
else 
   eegdir = [subdir '/eeg'];
end

if strcmp(session,'T1')
    sesdir = [eegdir '/ses-1'];
    if ~exist(sesdir)
        mkdir(sesdir);
    end
elseif strcmp(session,'T2')
    sesdir = [eegdir '/ses-2'];
    if ~exist(sesdir)
        mkdir(sesdir);
    end
elseif strcmp(session,'T3')
    sesdir = [eegdir '/ses-3'];
    if ~exist(sesdir)
        mkdir(sesdir);
    end
else
    fprintf('session directory error\n');
    sesdir = [eegdir '/error'];
end 

end

