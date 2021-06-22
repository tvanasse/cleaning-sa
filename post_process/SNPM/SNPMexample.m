netfile=[cd '/HydroCelGSN256v10.sfp'];
EEG.chanlocs = readlocs(netfile);
%insidegoodch = find([EEG.chanlocs(1:256).radius]<.57 & ~strcmp({EEG.chanlocs(1:256).labels},'Cz'));
insidegoodch = 1:256;
load([cd '/NeighborMatrix_256']);
neighbors(neighbors == 257) = NaN;
neighbors(257,:)=[];

%% Single or TFCE
E = 0.5;
H = 2;
alpha =.05;
comparison = 'pairedT';
tail = 'both'; % chose because you a priori believe post to have increased alpha
permutation_overide = 10000;  % 10000 permutations just to speed things up.


% straight up single threshold and TFCE test
[T,p,~,~]          = snpm_single_threshold_with_TFCE_NEW_FAST(data_x,data_y,neighbors,E,H,alpha,comparison,tail,permutation_overide);
sigch = find(p.corrected < .001); %sigch = find(p.corrected < alpha);
TopoplotSignificant(T.real_T,sigch,EEG.chanlocs,insidegoodch,'Single Threshold');

% cluster test using an abitrary critical T value = to uncorrected significance
degrees_of_freedom = number_of_subjects-2; % degrees of freedom for 19 paired subjects
%threshold          = 7; % very high threshold to show example of clustering normally would use line 34
threshold = tinv(1 - alpha/2,degrees_of_freedom); % if two sided divide alpha by 2;
[Clusters]         = snpm_cluster_NEW_FAST(data_x,data_y,threshold,neighbors,alpha,comparison,tail,permutation_overide);
sigclusters        = find([Clusters.p] < alpha);
sigch              = [Clusters(sigclusters).channels];
TopoplotSignificant(T.real_T,sigch,EEG.chanlocs,insidegoodch,'Clusters')



corr_x = [D_r2;ND_r2]; %
figure;topoplot(nanmean(corr_x),EEG.chanlocs,'maplimits','maxmin');title('ND_r2');colorbar
corr_y = repmat(DRS,1,256);
corr_x(:,outsidebad) = NaN;
corr_y(:,outsidebad) = NaN;
E = 0.5;
H = 2;
alpha =.05;
comparison = 'correlation';
tail = 'both'; % chose because you a priori believe post to have increased alpha
permutation_overide = 10000;  % 10000 permutations just to speed things up.

[corrT,corrp,~,~]   = snpm_single_threshold_with_TFCE_NEW_FAST(corr_x,corr_y,neighbors,E,H,alpha,comparison,tail,permutation_overide);
corrsigch           = find(corrp.corrected < .001); %sigch = find(p.corrected < alpha);
TopoplotSignificant(corrT.real_T,corrsigch,EEG.chanlocs,insidegoodch,'Single Threshold DRS correlation with delta');

% cluster test using an abitrary critical T value = to uncorrected significance
%degrees_of_freedom = 54; % degrees of freedom for 19 paired subjects
corrthreshold          = 0.5; % abitrary correlation threshold - could calculate a critical value based on degrees of freedom
[CorrClusters]         = snpm_cluster_NEW_FAST(corr_x,corr_y,corrthreshold,neighbors,alpha,comparison,tail,permutation_overide);
corrsigclusters        = find([CorrClusters.p] < alpha);
corrsigch              = [CorrClusters(corrsigclusters).channels];
TopoplotSignificant(corrT.real_T,corrsigch,EEG.chanlocs,insidegoodch,'Clusters DRS correlation with delta');

                        
   