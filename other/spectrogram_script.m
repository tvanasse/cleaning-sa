
EEG_all = pop_loadset('nrem_awakening_eeg_hp_trim_merged_nobadch_interp_avgref_ica_subcomps.set');
%EEG_all = pop_loadset('nrem_awakening_eeg_hp_trim_merged_nobadch_interp_avgref_ica.set');
addpath('../../../scripts/other/');

sesdir = pwd;

%nrem_awakenings = EEG_all.pnts/(EEG_all.srate*60*5 - EEG_all.srate*2 + 1)
save('nrem_index.mat','nrem_index');
load('nrem_index.mat');


for awak = 1:length(nrem_index)
    
    %get entire five minutes
    start_sample = (awak-1)*(EEG_all.srate*60*5 - EEG_all.srate*2 + 1);
    end_sample = awak*(EEG_all.srate*60*5 - EEG_all.srate*2 + 1); %each extraction is 4 min. 58seconds

    EEG = pop_select(EEG_all, 'point', [start_sample, end_sample]);

    EEG = pop_eegfiltnew(EEG, [], 50, [], 0, [], 0); %50 Hz low-pass filter

    subdir = ['awakening-' num2str(nrem_index(awak)) '-spectrograms'];
    if ~exist(subdir, 'dir')
       mkdir(subdir);
    end
    
    chan_spects = zeros(40,20,EEG.nbchan);
    for chan = 1:10
        dat = EEG.data(chan,EEG.pnts - EEG.srate*20:EEG.pnts);
        freq = [1:1:40];
        
        % 5-second epochs, 90% overlap
        [s,w,t] = spectrogram(dat,5*EEG.srate, 0.9*5*EEG.srate,freq,EEG.srate,'yaxis'); %channel vector, window sample size (in samples not seconds), sliding window size (i.e., % overlap), freq range, sampling rate
        chan_spects(:,:,chan) = abs(s);
        
        spectrogram(dat,1000,0,freq,500,'yaxis')
        %spectrogram(dat,500,250,[1:0.5:40],500,'yaxis')
        saveas(gcf, [subdir '/chan-' num2str(chan) '.png'], 'png')
        
        close all;
        
        % filter image with neighboring pixels
%         x = abs(s);
%         Iblur = imfilter(x,ones(3)/9)
%         imagesc(Iblur)
%         ax = gca;
%         ax.YDir = 'normal'
%         colorbar;
        
    end 
    
  
    
    save([subdir '/chan_spects.mat'],'chan_spects');

    EEG = pop_importdata('dataformat','array','data',EEG.data,...
        'srate',500,'xmin',0,'nbchan',EEG_all.nbchan, 'chanlocs', EEG_all.chanlocs);
    
    EEG = pop_saveset(EEG, 'filename', sprintf('awakening-%d-cleaned_nrem',nrem_index(awak)),'filepath',sesdir);
    
    psd_epoch_length = 6; %seconds
    upperHzlimit = 40; %Hz
    averef = 1;
    [psd,Hzbins] = psddata(EEG.data,EEG.srate,psd_epoch_length,upperHzlimit,averef);

    delta_idx = find(Hzbins > 1 & Hzbins <= 4);
    theta_idx = find(Hzbins > 4 & Hzbins <= 8);
    alpha_idx = find(Hzbins > 8 & Hzbins <= 12);
    beta_idx = find(Hzbins >= 12.5 & Hzbins <= 30);
    gamma_idx = find(Hzbins > 25 & Hzbins <= 40);

    freq_bans = {delta_idx,theta_idx,alpha_idx,beta_idx,gamma_idx};

    raw_topo = zeros(EEG.nbchan,5);
    ztopo = zeros(EEG.nbchan,5);

    for freq_l = 1:length(freq_bans)

        % psd -> channels x frequency_bins x epochs
        %Average across frequency bins between designations
        temp_topo = squeeze(mean(psd(:,freq_bans{freq_l},:),2));

        % average across all six-second epochs (no overlap)
        raw_topo(:,freq_l) = squeeze(mean(temp_topo(:,:),2));

        % z-score 
        ztopo(:,freq_l) = zscore(raw_topo(:,freq_l));

    end
    
    figure('Renderer', 'painters', 'Position', [100 100 1500 1500],'Name', sprintf('Awakening-%d-cleaned_nrem, 5 min.',nrem_index(awak)))
    hold on;
    ax1 = subplot(4,5,1); topoplot(raw_topo(:,1), EEG.chanlocs,'maplimits',[min(raw_topo(:,1)),max(raw_topo(:,1))],'electrodes','on','style','map'); title('Delta'); colorbar; colormap(ax1,jet)
    ax2 = subplot(4,5,2); topoplot(raw_topo(:,2), EEG.chanlocs,'maplimits',[min(raw_topo(:,2)),max(raw_topo(:,2))],'electrodes','on','style','map'); title('Theta '); colorbar; colormap(ax2,jet)
    ax3 = subplot(4,5,3); topoplot(raw_topo(:,3), EEG.chanlocs,'maplimits',[min(raw_topo(:,3)),max(raw_topo(:,3))],'electrodes','on','style','map'); title('Alpha'); colorbar; colormap(ax3,jet)
    ax5 = subplot(4,5,4); topoplot(raw_topo(:,4), EEG.chanlocs,'maplimits',[min(raw_topo(:,4)),max(raw_topo(:,4))],'electrodes','on','style','map'); title('Beta'); colorbar; colormap(ax5,jet)
    ax8 = subplot(4,5,5); topoplot(raw_topo(:,5), EEG.chanlocs,'maplimits',[min(raw_topo(:,5)),max(raw_topo(:,5))],'electrodes','on','style','map'); title('Gamma'); colorbar; colormap(ax8,jet)

    
    saveas(gcf, sprintf('awakening-%d-cleaned_nrem.tif',nrem_index(awak)), 'tif');
    %saveas(gcf, sprintf('awakening-%d-UN-cleaned_nrem.tif',nrem_index(awak)), 'tif');
    
    close all


    
end
