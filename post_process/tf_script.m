%% create multi-taper time x frequency plots/data


T = readtable('nrem_dataframe.csv');

tf_output = [];

fprintf('Total Files: %d\n',length(T.PATH));

for file = 1:length(T.PATH)
    output = [];
    
    fprintf('Iteration: %d\n',file);
    
    EEG = pop_loadset(char(T.PATH(file)));
    EEG = pop_select(EEG, 'point', [EEG.data - EEG.srate*120:length(EEG.data)-1]); % get 120 seconds before awakening

%     channel_index  = {'1','21','101','59','183','125','138'};
%     [~,~,chi]       = intersect(channel_labels,{EEG.chanlocs.labels});
    movingwin       = [5 .25];%%[2.5 .05];
    params.tapers   = [3 5];%[3 5];
    params.Fs       = EEG.srate;
    params.fpass    = [1 40];
    params.trialave = 0;
%     data            = EEG.data(sort(chi),:)';%EEG.data(sort(chi),:)';
    data            = EEG.data(:,:)';%EEG.data(sort(chi),:)';
    [S.mtspec,S.t,S.f]     = mtspecgramc(data,movingwin,params);

    % S.t is time points (seconds)
    % S.f is frequencies
    % S.mtspec is matrix with dim. len(S.t) x len(S.f) x channels
    minutes = S.t/(60);
    h = fspecial('disk',6);
    
    figure
    hmap = heatmap(S.t/60,S.f,log(imfilter(squeeze(S.mtspec(:,:,3)),h)'),...
                    'Colormap', parula,'GridVisible','off');
  

%     matrix = log(imfilter(squeeze(S.mtspec(:,:,3)),h)');
%     hmap = heatmap(downsample(S.t/60,5),...
%                     downsample(S.f,5), ...
%                     transpose(downsample(transpose(downsample(matrix,5)),5)),...
%                     'Colormap', parula,'GridVisible','off');
                
    hmap.YDisplayData = flip(hmap.YDisplayData);
    caxis([-6 4])

    % only show multples of ten
    idx = round(mod(S.t,10),2)==0; % index of datetime values to show as x tick
    hmap.XDisplayLabels(~idx) = {''}; % replace rejected tick labels with empties

    % only show multples of ten
    idx = round(mod(S.f,10),1)==0; % index of datetime values to show as x tick
    hmap.YDisplayLabels(~idx) = {''}; % replace rejected tick labels with empties

    xlabel('Time (Minutes)');
    ylabel('Frequency');
    set(gca,'fontsize',12);
    
    
    output.file = char(T.PATH(file));
    output.tf = transpose(downsample(transpose(downsample(matrix,5)),5));
    output.mtparams = params;
    output.length = length(EEG.data);
    tf_output = [output tf_output];
    
    save('mtspecrogram.mat','tf_output','-v7.3');
    
end
