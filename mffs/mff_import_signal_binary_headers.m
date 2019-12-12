function block_headers = mff_import_signal_binary_headers(file)
fid = fopen(file, 'r', 'l');


num_blocks   = 1;
block_headers       = struct([]);

%  fwrite(fidwrite, block_version, 'uint32');

while ~feof(fid)
   
    block_headers(num_blocks).block_version = fread(fid, 1, 'uint32'); % 0 
    % loop only done if there is a actual header
    
    if feof(fid)
        block_headers(num_blocks) = [];
        break
    end
    
    if block_headers(num_blocks).block_version > 0
        
        block_headers(num_blocks).header_size     = fread(fid, 1, 'uint32'); % 1
        % fwrite(fidwrite, block_header_size, 'uint32');  
        block_headers(num_blocks).data_size       = fread(fid, 1, 'uint32');
        % fwrite(fidwrite, block_data_size, 'uint32');
        block_headers(num_blocks).num_channels   = fread(fid, 1, 'uint32');
        % fwrite(fidwrite, block_num_channels, 'uint32');
        block_headers(num_blocks).channel_offsets = fread(fid, [1, block_headers(num_blocks).num_channels], 'uint32');
        % fwrite(fidwrite, block_channel_offsets, 'uint32');
        
        for channel_num = 1:block_headers(num_blocks).num_channels
            block_headers(num_blocks).sampling_depth(channel_num) = fread(fid, 1, 'uint8');
            %  fwrite(fidwrite, sampling_depth, 'uint8');
            
            if block_headers(num_blocks).sampling_depth(channel_num) ~= 32
                error('unsupported %d-bit sampling depth for channel %d in block %d', sample_size, channel_num, num_blocks);
            end
            
            block_headers(num_blocks).sampling_rate(channel_num) = fread(fid, 1, 'ubit24');
            % fwrite(fidwrite, sampling_rate, 'ubit24');
            
        end
        
        block_headers(num_blocks).optional_header_size = fread(fid, 1, 'uint32');
        % fwrite(fidwrite, optional_header_size, 'uint32');
        block_headers(num_blocks).optional_header = fread(fid, block_headers(num_blocks).optional_header_size);
        % fwrite(fidwrite, optional_header);
        
        fseek(fid, block_headers(num_blocks).data_size, 'cof');
        current_block_size = block_headers(num_blocks).data_size;
    else
        fseek(fid, current_block_size, 'cof');
    end
     num_blocks = num_blocks + 1;  
end

fclose(fid);

end