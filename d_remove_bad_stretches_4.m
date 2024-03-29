% ------------------------------------------------------------------------
% Plot ICA activations and remove bad stretches

% Author: Thomas Vanasse
% Center for Sleep and Consciousness, University of Wisconsin - Madison
% ------------------------------------------------------------------------

addpath('other/');
addpath('functions/');

eeglab;
close;

%% input filenames
inputlist = get_ses_dirs();

%% loop through input file list
for mff_input_file = 1:length(inputlist)
    
    sesdir = char(inputlist(mff_input_file));
    
    addpath(genpath(sesdir));
    
    load([sesdir '/nrem_index'])
    
    EEG_all = pop_loadset([sesdir '/nrem_awakening_eeg_hp_trim_merged_nobadch_ica.set']);
    
    if isempty(EEG_all.icaact) 
        EEG_all.icaact = (EEG_all.icaweights*EEG_all.icasphere)*EEG_all.data(EEG_all.icachansind,:);
    end
    
    for awak = 1:length(nrem_index)
        % print awakening
        fprintf('\n AWAKENING-%d, bad stretch removal\n\n',nrem_index(awak));

        % get entire five minutes
        start_sample = (awak-1)*(EEG_all.srate*60*5 - EEG_all.srate*2 + 1);
        fprintf('\n Start Sample (seconds): %d',start_sample/EEG_all.srate);
        
        end_sample = awak*(EEG_all.srate*60*5 - EEG_all.srate*2 + 1); %each extraction is 4 min. 58seconds
        fprintf('\n End Sample (seconds): %d',end_sample/EEG_all.srate);

        EEG = pop_select(EEG_all, 'point', [start_sample, end_sample]);
        
        EEG = pop_importdata('dataformat','array','data',EEG.data,...
        'srate',500,'xmin',0,'nbchan',EEG_all.nbchan, 'chanlocs', EEG_all.chanlocs);
    
        % get ica weights
        EEG.icaweights = EEG_all.icaweights;
        if awak == 1
            EEG.icaact = EEG_all.icaact(:, 1:end_sample);
        else 
            EEG.icaact = EEG_all.icaact(:, start_sample:end_sample);
        end
    
        EEG.icasphere = EEG_all.icasphere;

        eegplot(EEG.icaact,...
        'srate', EEG.srate,...
        'eloc_file', EEG.chanlocs,...
        'spacing', 15,...
        'winlength', 30,...
        'dispchans', 20,...
        'command', 'disp(''No data rejected. Use pop_select for this.'')',...
        'butlabel', 'MARK', ...
        'events', EEG.event);

        waitfor(gcf) %wait for epochs to be clean

        eegh %history function

        close all % close eeglab
        
        if size(TMPREJ,1) > 0
            % data storing the removed sections/bad channels
            EEG.badsections = TMPREJ(:,1:2);
            badsections = EEG.badsections;
            save([sesdir '/awakening-' num2str(nrem_index(awak)) '-badsections.mat'], 'badsections');

            % calculate percentage of removed sections
            sum_rmv = 0;
            for rmv = 1:size(TMPREJ,1)
                sum_rmv = sum_rmv + (TMPREJ(rmv,2)-TMPREJ(rmv,1));
                fprintf(num2str(rmv));
            end

            fid = fopen([sesdir '/log.txt'], 'at');
            fprintf(fid, 'Percentage of data removed, awakening-%d: %.2f\n', nrem_index(awak), sum_rmv/EEG.pnts*100);
            fclose(fid);
        
            % remove bad streches, but keep badsections data in 'badsections' field
            EEG = epi_log(@pop_select, EEG, 'nopoint', TMPREJ(:,1:2));
        
        end 

        EEG = pop_saveset(EEG, 'filename', sprintf('awakening-%d-nrem-posICA1-cleaned',nrem_index(awak)),'filepath',sesdir);
    
    end
    
end

