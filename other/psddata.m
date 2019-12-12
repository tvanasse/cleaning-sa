% 05/01/15 modified to not always take the average reference
function [psd,Hzbins] = psddata(data,sampling_rate,psd_epoch_length_in_secs,upperHzlimit,averef)
tic
if nargin <4 || isempty(upperHzlimit)
    upperHzlimit = sampling_rate/2;
end
psd_window_length = sampling_rate * psd_epoch_length_in_secs;

num_channels    = size(data,1);
num_samples     = size(data,2);
num_epochs      = floor(num_samples/psd_window_length);

start_samps     = 1:psd_window_length:num_samples;
end_samps       = psd_window_length:psd_window_length:num_samples;
Hzbins          = 0:1/psd_epoch_length_in_secs:upperHzlimit-1/psd_epoch_length_in_secs;

data_ave    = data(:,1:end)- repmat(nanmean(data(:,1:end)),num_channels,1);
if nargin == 5 && averef == 1 % data is already average referenced so don't do it
    data_ave    = data;
end
clear data;

psd  = NaN(num_channels,length(Hzbins),num_epochs);
for ep = 1:num_epochs
    display(ep)
        temp_psd    = pwelch(data_ave(:,start_samps(ep):end_samps(ep))',[],[],psd_window_length,sampling_rate);
        psd(:,:,ep) = temp_psd(1:length(Hzbins),:)';
        clear temp_psd;
end;
clear ep
% average across epochs
disp(['... elapsed time for computing spectral density ',num2str(toc/60),' minutes']);