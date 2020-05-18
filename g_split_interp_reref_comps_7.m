% ------------------------------------------------------------------------
% Get ICA cleaned data, split back to separate awakenings, then interpolate
% bad channels, and average reference data.

% Author: Thomas Vanasse
% Center for Sleep and Consciousness, University of Wisconsin - Madison
% ------------------------------------------------------------------------

addpath('other/');
addpath(genpath('functions/'));
eeglab;
close;

% start_dir = 0;

%% input filenames
more_files = 'Yes';
first_iter = 0;
while strcmp(more_files, 'No') ~= 1
batch_folder = uigetfile_n_dir();
if first_iter == 0;
    inputlist = uigetfile_n_dir(batch_folder);
else
    inputlist = [inputlist uigetfile_n_dir(batch_folder)];
end

more_files = questdlg('Add more subjects from different batch folder?', ...
	'??', ...
	'Yes','No','No');
first_iter = first_iter + 1;
end

% process earlier batch
% start_dir = '/Volumes/data-2/NCCAM3/SA/wDreamReport/aligned';
% inputlist = uigetfile_n_dir(start_dir);

%% loop through input file list

for mff_input_file = 1:length(inputlist)
    
    TABLE = readtable('NCCAM3_06_SADreamReports_10-20-18.csv','ReadVariableNames',true); %read experimenter input file
    filename = char(inputlist(mff_input_file)); %mff raw data filename
    
    k = strfind(filename, 'NCCAM_'); % get start index for filename (from full path)
    subid = filename(k+6:k+9);
    session = filename(k+11:k+12);
    
    %% create directory if it does not already exist
    current_dir = pwd;
    subdir = [current_dir '/../sub-' subid];
    
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
    
    %% plot all awakenings and sleep scoring data
    if isfile([sesdir '/all_awakenings.png'])
        
            [X,cmap] = imread([sesdir '/all_awakenings.png']);
            imshow(X,cmap);
            export_fig([sesdir '/' 'qa_figs.tiff'],'-append')
            close all;
    end
    
    %% plot good/bad components
    EEG = pop_loadset([sesdir '/nrem_merged_ica2.set']);
    
    % average reference for ICLabel
    EEG_avgref_beforeica2 = epi_log(@pop_reref, EEG, []);
    
    EEG_avgref_beforeica2 = iclabel(EEG_avgref_beforeica2); % prior to computing features, each dataset was converted to a common average reference
    
    % load badcomps
    badcomps = load([sesdir '/ic_artifacts.mat']);
    
    good_comps = [];
    for x = 1:size(EEG_avgref_beforeica2.icaweights,1)
        if isempty(find(badcomps.ic_artifacts_all==x))
            good_comps = [good_comps x];
        end      
    end
    
    pop_viewprops(EEG_avgref_beforeica2, 0, good_comps);
    
    figHandles = findobj('Type', 'figure');
    for x = size(figHandles,1):-1:1
        export_fig([sesdir '/' 'qa_figs.tiff'],'-append',figHandles(x))
    end
    
    close all;
    
    pop_viewprops(EEG_avgref_beforeica2, 0, badcomps.ic_artifacts_all);

    figHandles = findobj('Type', 'figure');
    for x = size(figHandles,1):-1:1
        export_fig([sesdir '/' 'qa_figs.tiff'],'-append',figHandles(x))
    end
    
    close all;
    
    [X,cmap] = imread([sesdir '/nrem_badchannels.tif']);
    imshow(X,cmap);
    export_fig([sesdir '/' 'qa_figs.tiff'],'-append')
    close all;
    
    %% continue with script
    load([sesdir '/nrem_index'])
    
    EEG_all = pop_loadset([sesdir '/nrem_merged_ica2_subcomps.set']);
    
    % remove additionally identified bad channels
%     EEG = epi_log(@pop_select, EEG, 'nochannel', EEG.badchannels);
    
    %50 Hz low-pass filter
    EEG_all = pop_eegfiltnew(EEG_all, [], 50, [], 0, [], 0); 
    
    load([sesdir '/cleaned_lengths.mat']);
    cumul_sample = 1; %keep track of sample 
    
for awak = 1:length(nrem_index)
    
    %get entire five minutes
    start_sample = cumul_sample;
    end_sample = cumul_sample + cleaned_lengths(awak) - 1; 
    fprintf(['start: ' num2str(start_sample) ' end: ' num2str(end_sample) '\n'])
    cumul_sample = cumul_sample + cleaned_lengths(awak);

    %% Plot topo plots before ica
    EEG = pop_select(EEG_avgref_beforeica2, 'point', [start_sample, end_sample]);

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
    spindle_idx = find(Hzbins > 12 & Hzbins <= 15);
    beta_idx = find(Hzbins >= 12.5 & Hzbins <= 30);
    gamma_idx = find(Hzbins > 25 & Hzbins <= 40);

    freq_bans = {delta_idx,theta_idx,alpha_idx,spindle_idx,beta_idx,gamma_idx};

    raw_topo = zeros(EEG_avgref.nbchan,6);
    ztopo = zeros(EEG_avgref.nbchan,6);

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
    ax1 = subplot(1,6,1); topoplot(raw_topo(:,1), EEG_avgref.chanlocs,'maplimits',[min(raw_topo(:,1)),max(raw_topo(:,1))],'electrodes','on','style','map'); title(sprintf('Delta-A%d-BICA', nrem_index(awak))); colorbar; colormap(ax1,jet)
    ax2 = subplot(1,6,2); topoplot(raw_topo(:,2), EEG_avgref.chanlocs,'maplimits',[min(raw_topo(:,2)),max(raw_topo(:,2))],'electrodes','on','style','map'); title(sprintf('Theta-A%d-BICA', nrem_index(awak))); colorbar; colormap(ax2,jet)
    ax3 = subplot(1,6,3); topoplot(raw_topo(:,3), EEG_avgref.chanlocs,'maplimits',[min(raw_topo(:,3)),max(raw_topo(:,3))],'electrodes','on','style','map'); title(sprintf('Alpha-A%d-BICA', nrem_index(awak))); colorbar; colormap(ax3,jet)
    ax4 = subplot(1,6,4); topoplot(raw_topo(:,4), EEG_avgref.chanlocs,'maplimits',[min(raw_topo(:,4)),max(raw_topo(:,4))],'electrodes','on','style','map'); title(sprintf('Spindles-A%d-BICA', nrem_index(awak))); colorbar; colormap(ax3,jet)
    ax5 = subplot(1,6,5); topoplot(raw_topo(:,5), EEG_avgref.chanlocs,'maplimits',[min(raw_topo(:,5)),max(raw_topo(:,5))],'electrodes','on','style','map'); title(sprintf('Beta-A%d-BICA', nrem_index(awak))); colorbar; colormap(ax5,jet)
    ax8 = subplot(1,6,6); topoplot(raw_topo(:,6), EEG_avgref.chanlocs,'maplimits',[min(raw_topo(:,6)),max(raw_topo(:,6))],'electrodes','on','style','map'); title(sprintf('Gamma-A%d-BICA', nrem_index(awak))); colorbar; colormap(ax8,jet)

    saveas(gcf, [sesdir '/' sprintf('awakening-%d-cleaned2_nrem.tif', nrem_index(awak))], 'tif');
    export_fig([sesdir '/' 'qa_figs.tiff'],'-append')
    close all
    
    %% Plot topo plots after ica
    EEG = pop_select(EEG_all, 'point', [start_sample, end_sample]);

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
    spindle_idx = find(Hzbins > 12 & Hzbins <= 15);
    beta_idx = find(Hzbins >= 12.5 & Hzbins <= 30);
    gamma_idx = find(Hzbins > 25 & Hzbins <= 40);

    freq_bans = {delta_idx,theta_idx,alpha_idx,spindle_idx,beta_idx,gamma_idx};

    raw_topo = zeros(EEG_avgref.nbchan,6);
    ztopo = zeros(EEG_avgref.nbchan,6);

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
    ax1 = subplot(1,6,1); topoplot(raw_topo(:,1), EEG_avgref.chanlocs,'maplimits',[min(raw_topo(:,1)),max(raw_topo(:,1))],'electrodes','on','style','map'); title(sprintf('Delta-A%d', nrem_index(awak))); colorbar; colormap(ax1,jet)
    ax2 = subplot(1,6,2); topoplot(raw_topo(:,2), EEG_avgref.chanlocs,'maplimits',[min(raw_topo(:,2)),max(raw_topo(:,2))],'electrodes','on','style','map'); title(sprintf('Theta-A%d', nrem_index(awak))); colorbar; colormap(ax2,jet)
    ax3 = subplot(1,6,3); topoplot(raw_topo(:,3), EEG_avgref.chanlocs,'maplimits',[min(raw_topo(:,3)),max(raw_topo(:,3))],'electrodes','on','style','map'); title(sprintf('Alpha-A%d', nrem_index(awak))); colorbar; colormap(ax3,jet)
    ax4 = subplot(1,6,4); topoplot(raw_topo(:,4), EEG_avgref.chanlocs,'maplimits',[min(raw_topo(:,4)),max(raw_topo(:,4))],'electrodes','on','style','map'); title(sprintf('Spindles-A%d', nrem_index(awak))); colorbar; colormap(ax3,jet)
    ax5 = subplot(1,6,5); topoplot(raw_topo(:,5), EEG_avgref.chanlocs,'maplimits',[min(raw_topo(:,5)),max(raw_topo(:,5))],'electrodes','on','style','map'); title(sprintf('Beta-A%d', nrem_index(awak))); colorbar; colormap(ax5,jet)
    ax8 = subplot(1,6,6); topoplot(raw_topo(:,6), EEG_avgref.chanlocs,'maplimits',[min(raw_topo(:,6)),max(raw_topo(:,6))],'electrodes','on','style','map'); title(sprintf('Gamma-A%d', nrem_index(awak))); colorbar; colormap(ax8,jet)

    saveas(gcf, [sesdir '/' sprintf('awakening-%d-cleaned2_nrem.tif', nrem_index(awak))], 'tif');
    export_fig([sesdir '/' 'qa_figs.tiff'],'-append')
    close all
    % Plot power spectra of awakening:
    mat_path = sesdir;
    mat_name = sprintf('awakening-%d-cleaned2_nrem.mat', nrem_index(awak));
    options = struct(...
        'save_file',        0            ,... 
        'save_name',        mat_name     ,...
        'save_path',        mat_path     , ...
        'epoch_length',     6            ,...  % window length in seconds
        'freq_limit',       240          ,...  % number of bands
        'ylimitmax',        1            , ... 
        'fft_bands',        1)           ;     

    % Calculate the FFT
    % Stores output variables and input options in specified file.
    [fft_all, freq_range] = csc_FFT_last(EEG_avgref, options);
    
    % "plot the spectra"
    Spectra_Plotter(fft_all,1:size(fft_all,3), options);
    export_fig([sesdir '/' 'qa_figs.tiff'],'-append')
    close all;
    
    [X,cmap] = imread([sesdir sprintf('/awakening-%d-scoring.png', nrem_index(awak))]);
    imshow(X,cmap);
    export_fig([sesdir '/' 'qa_figs.tiff'],'-append')
    close all;
    
  
end
end
