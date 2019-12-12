%Manually remove bad channels, interpolate, and average reference
CHANNEL_LOCATION_FILE = 'channel_location_file/HydroCelGSN256v10.sfp';

selpath = uigetdir;
set_files = dir([selpath '/*.set']);
set_files = size(set_files,1)/2;

for i = 1:set_files
EEG = pop_loadset([selpath '/awakening-' num2str(i) '_eeg_hp.set']);

EEG = csc_eeg_plotter(EEG);

%fig = gcf; %Some Figure
waitfor(EEG);

EEG.badchannels = sort(EEG.csc_hidden_channels);

figure; topoplot([],EEG.chanlocs([EEG.badchannels]), 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo); % which channels removed
saveas(gcf, [selpath '/' 'awakening-' num2str(i) '_badchan.tif'], 'tif');
close all;

% Remove the bad channels, saving original channel info for later
EEG = epi_log(@pop_select, EEG, 'nochannel', EEG.badchannels);

EEG.urchanlocs = readlocs(CHANNEL_LOCATION_FILE);

%EEG.urchanlocs = EEG.chanlocs; 
% Interpolate.
EEG = epi_log(@eeg_interp, EEG, EEG.urchanlocs); 

% Average reference
EEG = epi_log(@pop_reref, EEG, []);

% Save to disk.
% Ben+Anna deleted the filtered_file, but I have opted to leave it on disk.
nobadch_interpol_avg_file = [EEG.filename(1:end-4) '_nobadch_interpol_avg'];
EEG = epi_log(@pop_saveset, EEG, 'filename', nobadch_interpol_avg_file,'filepath',selpath, 'check', 'on');

end
