netfile='../channel_location_file/HydroCelGSN256v10.sfp';
addpath('SNPM');
EEG.chanlocs = readlocs(netfile);

% insidegoodch = 173 channels
insidegoodch = find([EEG.chanlocs(1:256).radius]<.57 & ~strcmp({EEG.chanlocs(1:256).labels},'Cz'))
% topoplot([],EEG.chanlocs(insidegoodch))


% insidegoodch = 1:185;
load(['SNPM/NeighborMatrix_256']);

neighbors(neighbors == 257) = NaN;
neighbors(257,:)=[];

% inside185ch = 185 channels
load('../channel_location_file/inside185ch.mat');
topoplot([],EEG.chanlocs(inside185ch))
EEG.chanlocs(inside185ch).labels

insidegoodch = inside185ch

% overlap, returns position in first argument, inside185ch
% [val,pos]=intersect(inside185ch,insidegoodch)

%% Load PSD Data
psd_dir = '/data/tvanasse/nccam3/data/post_process_data/psds/';

for i={'delta','theta','alpha','beta','gamma'}
   
disp(i{1})

data_x = load([psd_dir 'reportandnoreport/report_freqnorm_60s_idavg.mat']);
data_x = eval(['data_x.' i{1}]);
%data_x = data_x(:,pos)
size(data_x)

data_y = load([psd_dir 'reportandnoreport/noreport_freqnorm_60s_idavg.mat']);
data_y = eval(['data_y.' i{1}]);

%data_y = data_y(:,pos)
size(data_y)

%% Load Spectral Exponents
% specexps_dir = '/data/tvanasse/nccam3/data/post_process_data/specexps/';
% data_x = load([specexps_dir 'reportorsomethingandnoreport/reportsomething_SPECEXP_1-40_ALL_idavg_60s.mat']);
% data_x = data_x.sub_chan;
% size(data_x)
% 
% data_y = load([specexps_dir 'reportorsomethingandnoreport/noreport_SPECEXP_1-40_ALL_idavg_60s.mat']);
% data_y = data_y.sub_chan;
% size(data_y)

%% Fill data outside 185's with NaN
% data_x(:,setdiff(1:256,inside185ch))=NaN;
% data_y(:,setdiff(1:256,inside185ch))=NaN;

data_x_256 = NaN(size(data_x,1),256);
data_x_256(:,inside185ch) = data_x;

data_y_256 = NaN(size(data_y,1),256);
data_y_256(:,inside185ch) = data_y;

size(data_x_256)
size(data_y_256)

%%
E = 0.5;
H = 2;
alpha =.05;

% addpath('/Volumes/apps/linux/R2020b/toolbox/matlab/strfun')
comparison = 'pairedT';
tail = 'right'; 
permutation_overide = 1000;  % 10000 permutations just to speed things up.

% % straight up single threshold and TFCE test
[T,p,~,~]          = snpm_single_threshold_with_TFCE_NEW_FAST(data_x_256,data_y_256,neighbors,E,H,alpha,comparison,tail,permutation_overide);
% sigch = find(p.corrected < .01); %sigch = find(p.corrected < alpha);
% sigch
% TopoplotSignificant(T.real_T,sigch,EEG.chanlocs,inside185ch,'Single Threshold');
% colormap default
% caxis([-2.5,2.5])
% 

% cluster test using an abitrary critical T value = to uncorrected significance
degrees_of_freedom = size(data_x,1) + size(data_y,1) - 2; 
%threshold          = 7; % very high threshold to show example of clustering normally would use line 34
threshold = tinv(1 - alpha/2,degrees_of_freedom); % if two sided divide alpha by 2;
threshold
fprintf('Threshold: %d \n', threshold);
[Clusters]         = snpm_cluster_NEW_FAST(data_x_256,data_y_256,threshold,neighbors,alpha,comparison,tail,permutation_overide);
sigclusters        = find([Clusters.p] < alpha);
sigch              = [Clusters(sigclusters).channels];
sigclusters
sigch
TopoplotSignificant(T.real_T,sigch,EEG.chanlocs,inside185ch,'Clusters')
colormap default
caxis([-2.5,2.5])



title(i{1})
end


% corr_x = [D_r2;ND_r2]; %
% figure;topoplot(nanmean(corr_x),EEG.chanlocs,'maplimits','maxmin');title('ND_r2');colorbar
% corr_y = repmat(DRS,1,256);
% corr_x(:,outsidebad) = NaN;
% corr_y(:,outsidebad) = NaN;
% E = 0.5;
% H = 2;
% alpha =.05;
% comparison = 'correlation';
% tail = 'both'; % chose because you a priori believe post to have increased alpha
% permutation_overide = 10000;  % 10000 permutations just to speed things up.
% 
% [corrT,corrp,~,~]   = snpm_single_threshold_with_TFCE_NEW_FAST(corr_x,corr_y,neighbors,E,H,alpha,comparison,tail,permutation_overide);
% corrsigch           = find(corrp.corrected < .001); %sigch = find(p.corrected < alpha);
% TopoplotSignificant(corrT.real_T,corrsigch,EEG.chanlocs,insidegoodch,'Single Threshold DRS correlation with delta');
% 
% % cluster test using an abitrary critical T value = to uncorrected significance
% %degrees_of_freedom = 54; % degrees of freedom for 19 paired subjects
% corrthreshold          = 0.5; % abitrary correlation threshold - could calculate a critical value based on degrees of freedom
% [CorrClusters]         = snpm_cluster_NEW_FAST(corr_x,corr_y,corrthreshold,neighbors,alpha,comparison,tail,permutation_overide);
% corrsigclusters        = find([CorrClusters.p] < alpha);
% corrsigch              = [CorrClusters(corrsigclusters).channels];
% TopoplotSignificant(corrT.real_T,corrsigch,EEG.chanlocs,insidegoodch,'Clusters DRS correlation with delta');
% 
                        
   