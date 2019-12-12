% make changes to reading breaks to account for bug
function [num_epoch_breaks, epoch_breaks] = mff_determine_epoch_breaks(num_epochs, epochs)
    num_epoch_breaks = 0;
    
    if num_epochs > 1
        num_epoch_breaks            = num_epochs - 1;
        for epoch_num = 1:num_epochs - 1
            epoch_breaks.onset(epoch_num)    = epochs.time_to(epoch_num);
            epoch_breaks.duration(epoch_num) = epochs.time_from(epoch_num + 1) - epochs.time_to(epoch_num);
        end
    else 
       epoch_breaks     = struct([]); 
    end
end
