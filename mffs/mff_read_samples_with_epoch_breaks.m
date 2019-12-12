function sampledata = mff_read_samples_with_epoch_breaks(meta_file,channels_index,sample_start, sample_end)
    meta_data = mff_import_meta_data(meta_file);
    
    num_channels           = meta_data.signal_binaries.num_channels;
     
    sample_index           = sample_start:sample_end;% - block_samples_starts(min_block_index) + 1;
    
    if isequal(channels_index,'all')
        channels_index = 1:num_channels;
    end
    
    sampledata=cell(length(channels_index),1);
    if max(channels_index > num_channels)
        error('channel indicies out of range of data')
    end
    i=1;
    for channel_num = channels_index
        channels = mff_import_signal_binary_with_epoch_breaks_allblocks(meta_data, channel_num); 
        sampledata{i}=channels.samples(sample_index);
        i=i+1;
    end
    sampledata = cell2mat(sampledata);
   
end