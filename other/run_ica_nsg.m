function [outputArg1] = run_ica_nsg(EEG,jobpath)
%EEG structure must have EEG.filename specified

eeglab;
close; 

% Runnning a job with default options by passing the path to the job and
% defining the script to execute inside the job's zip file

% File writing begin ---
fid = fopen( fullfile(jobpath,'icansg_job.m'), 'w');
fprintf(fid, 'eeglab;\n');
fprintf(fid, 'EEG = pop_loadset(''%s'');\n', EEG.filename);
fprintf(fid, 'EEG = pop_runica(EEG, ''%s'',''%s'');\n', 'icatype','runica');
fprintf(fid, 'pop_saveset(EEG, ''filename'', ''%s'');\n',[EEG.filename(1:end-4) '_ica.set']);
fclose(fid);
% File writing end ---

zip('nsg_output',jobpath)

path2zip = '/Volumes/data/NCCAM3/SA/wDreamReport/aligned/extraction_TJV/sub-2052/eeg/ses-3/nsg_output.zip';


nsg_jobstruct = pop_nsg('run',path2zip,'filename', 'icansg_job.m','jobid',EEG.filename,'runtime',3);

% Alternatively you can submit a job and provide extra options as 'runtime', 'jobid'.
% [currentjob1, alljobs] = pop_nsg('run',path2zip, 'jobid', 'runica_testing','runtime', 0.3 ); 

%% Checking job status periodically
% Check the status of the job periodically by calling the function nsg_recurspol, providing 
% as argument the NSG job structure(see example below), a job ID or a job URL

jobstructout = nsg_recurspoll('run_ica_test','pollinterval', 30);

% jobstructout = nsg_recurspoll(jobstruct,...
%                               'pollinterval', POLL_INTERVAL); % recursively poll for NSG job status 
%                                                               % and display latest results

%% Retrieving job results
% Retreive and download job results by providing an NSG job structure (see example below),
% a job ID or a job URL to pop_nsg
[currentjobout, alljobs] = pop_nsg('output',jobstruct); 

%% Deleting an NSG job
% Delete a job from the NSG record associated with the user NSG credential by
% providing an NSG job structure (see example below), a job ID or a job URL to pop_nsg

[jobdeleted, alljobs] = pop_nsg('delete','runica_testing');

outputArg1 = 1;

end

