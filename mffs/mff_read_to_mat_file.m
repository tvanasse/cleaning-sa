function matdata = mff_read_to_mat_file(meta_file)
meta_data = mff_import_meta_data(meta_file);

num_channels = meta_data.signal_binaries(1).num_channels;
num_samples  = meta_data.signal_binaries(1).channels.num_samples(1);

% warning('ONLY ONE CHANNEL');
% num_channels = 1;

matdata=NaN(num_channels,num_samples);

for channel_num = 1:num_channels
    
    if ~mod(channel_num, round(num_channels/10))
        fprintf('.')
        drawnow;
    end
    
    channels = mff_import_signal_binary_with_epoch_breaks_allblocks(meta_data, channel_num); 
    
    matdata(channel_num,:) = channels(1).samples;    
end

samplingRate = meta_data.signal_binaries(1).channels.sampling_rate(1); %#ok<NASGU>

end
