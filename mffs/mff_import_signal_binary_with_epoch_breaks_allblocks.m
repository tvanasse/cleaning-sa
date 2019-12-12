function channels = mff_import_signal_binary_with_epoch_breaks_allblocks(mff_meta_data, channels_index)
% write function to deal with multiple signal binaries
if isequal(channels_index, 'all')
    channels_index = 1:mff_meta_data.signal_binaries.num_channels;
end

num_channels = length(channels_index);
num_blocks = mff_meta_data.signal_binaries.num_blocks;%length(blocks_list);
blocks_list = 1:num_blocks;

id = fopen(mff_meta_data.signal_binaries.file, 'r', 'l');

block_offsets          = mff_meta_data.signal_binaries.blocks.offset(blocks_list, channels_index);
block_num_samples      = mff_meta_data.signal_binaries.blocks.num_samples(blocks_list, channels_index);
% channel_offsets        = [ones(1, num_channels); cumsum(block_num_samples(1:end-1, :)) + 1];
channel_sampling_rates = mff_meta_data.signal_binaries.channels.sampling_rate(channels_index);
channel_num_samples    = round((mff_meta_data.epochs.time_to(end)- mff_meta_data.epochs.time_from(1))*channel_sampling_rates);%sum(block_num_samples, 1);
calibrated_gains       = mff_meta_data.signal_binaries.calibrated_gains(blocks_list, channels_index);
calibrated_zeros       = mff_meta_data.signal_binaries.calibrated_zeros(blocks_list, channels_index);

channels = repmat(struct('num', NaN, 'sampling_rate', NaN, 'samples', []), num_channels, 1);

for channel_num = 1:num_channels
    channels(channel_num).num           = channels_index(channel_num);
    channels(channel_num).sampling_rate = channel_sampling_rates(channel_num);
    channels(channel_num).samples       = zeros(1, channel_num_samples(channel_num));
end

for epoch_num = 1:mff_meta_data.num_epochs
    
    sample_from = round(mff_meta_data.epochs.time_from(epoch_num)*channels(channel_num).sampling_rate+1);%channel_offsets(block_num, channel_num); % use mff_meta_data.time_from
     
    for block_num = mff_meta_data.epochs.block_from(epoch_num):...
                    mff_meta_data.epochs.block_to(epoch_num)
        
        for channel_num = 1:num_channels
            offset = block_offsets(block_num, channel_num);
            
            fseek(id, offset, 'bof');
            
            num_samples = block_num_samples(block_num, channel_num);
            samples     = fread(id, [1, num_samples], '*single');
            
            sample_to   = sample_from + num_samples - 1;
            
             channels(channel_num).samples(sample_from:sample_to) = ...
                 calibrated_gains(block_num, channel_num) * ...
                 (samples(:) - calibrated_zeros(block_num, channel_num));
%             
        end
        % update samples to reflect new block
        sample_from = sample_from + num_samples;
    end
end
fclose(id);
% break data appears to be inaccurate
for channel_num = 1:num_channels
    for break_num = 1:mff_meta_data.num_epoch_breaks
        num_break_samples  = round(mff_meta_data.epoch_breaks.duration(break_num)*channels(channel_num).sampling_rate);
        break_from   = mff_meta_data.epoch_breaks.onset(break_num)*channels(channel_num).sampling_rate+1;
        break_to     = break_from + num_break_samples - 1;
        breaksamples = linspace(channels(channel_num).samples(break_from-1),...
            channels(channel_num).samples(break_to+1),num_break_samples);
        channels(channel_num).samples(break_from:break_to) = breaksamples;
    end
end
end
