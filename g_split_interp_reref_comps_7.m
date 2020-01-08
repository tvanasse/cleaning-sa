% ------------------------------------------------------------------------
% Get ICA cleaned data, split back to separate awakenings, then interpolate
% bad channels, and average reference data.

% Author: Thomas Vanasse
% Center for Sleep and Consciousness, University of Wisconsin - Madison
% ------------------------------------------------------------------------

addpath('other/');
addpath('functions/');

inputlist = uigetfile_n_dir;

%% loop through input file list

for mff_input_file = 1:length(inputlist)
    
    TABLE = readtable('NCCAM3_06_SADreamReports_10-20-18.csv','ReadVariableNames',true); %read experimenter input file
    filename = char(inputlist(mff_input_file)); %mff raw data filename
    
    k = strfind(filename, 'NCCAM_'); % get start index for filename (from full path)
    subid = filename(k+6:k+9);
    session = filename(k+11:k+12);
    
    %% create directory if it does not already exist
    subdir = ['/Volumes/data/NCCAM3/SA/wDreamReport/aligned/extraction_TJV' '/sub-' subid];
    if ~exist(subdir, 'dir')
       mkdir(subdir);
       eegdir = [subdir '/eeg'];
       mkdir(eegdir);
    else 
       eegdir = [subdir '/eeg'];
    end
    
    if strcmp(session,'T1')
        sesdir = [eegdir '/ses-1'];
        timepoint=1;
        if ~exist(sesdir)
            mkdir(sesdir);
        end
    elseif strcmp(session,'T2')
        sesdir = [eegdir '/ses-2'];
        timepoint=2;
        if ~exist(sesdir)
            mkdir(sesdir);
        end
    elseif strcmp(session,'T3')
        sesdir = [eegdir '/ses-3'];
        timepoint=3;
        if ~exist(sesdir)
            mkdir(sesdir);
        end
    else
        fprintf('session directory error\n');
        sesdir = [eegdir '/error'];
    end 
    
    load([sesdir '/nrem_index'])
    
    EEG_all = pop_loadset([sesdir '/nrem_merged_ica2_subcomps.set']);
    
    load([sesdir '/cleaned_lengths.mat']);
    cumul_sample = 1; %keep track of sample 
    
for awak = 1:length(nrem_index)
    
        %get entire five minutes
        start_sample = cumul_sample;
        end_sample = cumul_sample + cleaned_lengths(awak) - 1; 
        cumul_sample = cumul_sample + cleaned_lengths(awak);

        EEG = pop_select(EEG_all, 'point', [start_sample, end_sample]);

%     subdir = ['awakening-' num2str(nrem_index(awak)) '-spectrograms'];
%     if ~exist(subdir, 'dir')
%        mkdir(subdir);
%     end
%     
%     chan_spects = zeros(40,20,EEG.nbchan);
%     for chan = 1:10
%         dat = EEG.data(chan,EEG.pnts - EEG.srate*20:EEG.pnts);
%         freq = [1:1:40];
%         
%         % 5-second epochs, 90% overlap
%         [s,w,t] = spectrogram(dat,5*EEG.srate, 0.9*5*EEG.srate,freq,EEG.srate,'yaxis'); %channel vector, window sample size (in samples not seconds), sliding window size (i.e., % overlap), freq range, sampling rate
%         chan_spects(:,:,chan) = abs(s);
%         
%         spectrogram(dat,1000,0,freq,500,'yaxis')
%         %spectrogram(dat,500,250,[1:0.5:40],500,'yaxis')
%         saveas(gcf, [subdir '/chan-' num2str(chan) '.png'], 'png')
%         
%         close all;
%         
%         % filter image with neighboring pixels
% %         x = abs(s);
% %         Iblur = imfilter(x,ones(3)/9)
% %         imagesc(Iblur)
% %         ax = gca;
% %         ax.YDir = 'normal'
% %         colorbar;
%         
%     end 
       
%     save([subdir '/chan_spects.mat'],'chan_spects');

    EEG = pop_importdata('dataformat','array','data',EEG.data,...
        'srate',500,'xmin',0,'nbchan',EEG_all.nbchan, 'chanlocs', EEG_all.chanlocs);

    % load bad sections
    
    if isfile([sesdir '/awakening-' num2str(nrem_index(awak)) '-badsections.mat'])
        load([sesdir '/awakening-' num2str(nrem_index(awak)) '-badsections.mat']);
        EEG.badsections = badsections;
    else
        EEG.badsections = []; % no bad sections
    end
    
    
    
    % identify and interpoloate bad channels
    load([sesdir '/chanlocs_185.mat'])
    EEG.urchanlocs = origEEGchanlocs;   
    EEG = epi_log(@eeg_interp, EEG, EEG.urchanlocs); 
    
    %50 Hz low-pass filter
    EEG = pop_eegfiltnew(EEG, [], 50, [], 0, [], 0); 
    
    % average reference
    EEG = epi_log(@pop_reref, EEG, []);

    EEG = pop_saveset(EEG, 'filename', sprintf('awakening-%d-cleaned2_nrem',nrem_index(awak)),'filepath',sesdir);
    
    %Calculate point spectral density
    % average reference
    EEG_avgref = EEG;
    psd_epoch_length = 6; %seconds
    upperHzlimit = 40; %Hz
    averef = 1;
    [psd,Hzbins] = psddata(EEG_avgref.data,EEG_avgref.srate,psd_epoch_length,upperHzlimit,averef);

    delta_idx = find(Hzbins > 1 & Hzbins <= 4);
    theta_idx = find(Hzbins > 4 & Hzbins <= 8);
    alpha_idx = find(Hzbins > 8 & Hzbins <= 12);
    beta_idx = find(Hzbins >= 12.5 & Hzbins <= 30);
    gamma_idx = find(Hzbins > 25 & Hzbins <= 40);

    freq_bans = {delta_idx,theta_idx,alpha_idx,beta_idx,gamma_idx};

    raw_topo = zeros(EEG_avgref.nbchan,5);
    ztopo = zeros(EEG_avgref.nbchan,5);

    for freq_l = 1:length(freq_bans)
        % psd -> channels x frequency_bins x epochs
        % average across frequency bins between designations
        temp_topo = squeeze(mean(psd(:,freq_bans{freq_l},:),2));

        % average across all six-second epochs (no overlap)
        raw_topo(:,freq_l) = squeeze(mean(temp_topo(:,:),2));

        % z-score 
        ztopo(:,freq_l) = zscore(raw_topo(:,freq_l));
    end
      
    figure('Renderer', 'painters', 'Position', [100 900 1500 200],'Name', sprintf('Awakening-%d-cleaned_nrem, 5 min.',nrem_index(awak)))
    hold on;
    ax1 = subplot(1,5,1); topoplot(raw_topo(:,1), EEG_avgref.chanlocs,'maplimits',[min(raw_topo(:,1)),max(raw_topo(:,1))],'electrodes','on','style','map'); title('Delta'); colorbar; colormap(ax1,jet)
    ax2 = subplot(1,5,2); topoplot(raw_topo(:,2), EEG_avgref.chanlocs,'maplimits',[min(raw_topo(:,2)),max(raw_topo(:,2))],'electrodes','on','style','map'); title('Theta '); colorbar; colormap(ax2,jet)
    ax3 = subplot(1,5,3); topoplot(raw_topo(:,3), EEG_avgref.chanlocs,'maplimits',[min(raw_topo(:,3)),max(raw_topo(:,3))],'electrodes','on','style','map'); title('Alpha'); colorbar; colormap(ax3,jet)
    ax5 = subplot(1,5,4); topoplot(raw_topo(:,4), EEG_avgref.chanlocs,'maplimits',[min(raw_topo(:,4)),max(raw_topo(:,4))],'electrodes','on','style','map'); title('Beta'); colorbar; colormap(ax5,jet)
    ax8 = subplot(1,5,5); topoplot(raw_topo(:,5), EEG_avgref.chanlocs,'maplimits',[min(raw_topo(:,5)),max(raw_topo(:,5))],'electrodes','on','style','map'); title('Gamma'); colorbar; colormap(ax8,jet)

    saveas(gcf, [sesdir '/' sprintf('awakening-%d-cleaned2_nrem.tif', nrem_index(awak))], 'tif');
    
    close all
  
end
end
