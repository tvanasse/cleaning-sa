%EEG = pop_loadset('/Volumes/NCCAM/NCCAM/NCCAM3/NCCAM3_Workspace_Backup/wDreamReport/aligned/extraction_TJV/cleaning-sa/post_process/../../sub-2001/eeg/ses-3/awakening-2-cleaned2_nrem.set');


% remove alternate path to findpeaks (used in fitPowerLaw)
rmpath('/Users/tononilab/Documents/MATLAB/eeglab14_1_2b/plugins/PrepPipeline0.55.3/utilities/chronux_2_modified/spectral_analysis/continuous/')


T = readtable('nrem_dataframe.csv');

output = [];

fprintf('Total Files: %d\n',length(T.PATH));
for file = 1:length(T.PATH)
    fprintf('Iteration: %d\n',file);
    
    EEG = pop_loadset(char(T.PATH(file)));
%     if strcmp(char(T.HOLD_OUT(file)),'False')
%         EEG = pop_loadset(char(T.PATH(file)));
%     else
%         continue
%     end

    
    for timerange = [1, 2] % minutes before awakening
        
        for freqranges = {[1,40], [1 20], [20 40]} % range of frequencies for slope
            
            chanavg = []; % average of spectral exponents
            obsavg = []; % average of observed psd
            
            for chan = 1:EEG.nbchan % 185 channels
                 %compute the PSD
                 epLen=2*EEG.srate; epShift=1*EEG.srate; numFFT=[];
                 [myPSD,frex]= pwelch( EEG.data(chan,size(EEG.data,2)-timerange*EEG.srate*60:size(EEG.data,2)), ...
                                       epLen, epShift,numFFT, EEG.srate); 

                 %frBand=[1 40];
                 frBand = freqranges{1};
                 frBins= dsearchn( frex, frBand(1) ):  dsearchn( frex, frBand(2));
                 XX= frex(frBins);
                 YY= myPSD(frBins);
                 robRegMeth= 'ols'; % method to perform linear regression. see >> help robustfit

                 %doPlot= 1; figure;
                 doPlot = 0;

                 thisCol= [0 0 1];
                 [intSlo, stat, Pows, Deviants,  stat0, intSlo0] = fitPowerLaw3steps(XX,YY, robRegMeth,  doPlot, thisCol);
                 spectralExponent= intSlo(2);
                 
                 chanavg = [spectralExponent chanavg];
                 obsavg = [Pows.obs ; obsavg];
                 
            end
                 
             %store output in structure
             data.path = T.PATH(file);
             data.chanlocs = EEG.chanlocs(chan);
%                  data.intSlo = intSlo;
%                  data.stat = stat;
%                  data.Pows = Pows;
             data.Deviants = Deviants;
%                  data.stat0 = stat0;
             data.intSlo0 = intSlo0;
             data.freqrange = frBand;
             data.timerange = timerange;
                 
             data.Pows = Pows;
             
             data.spectralexp_all = chanavg; % spectral exponent for each channel
             data.obs_freq_all = obsavg; % observed frequencies for each channel (matrix)
             
             data.meanspectralexp = mean(chanavg);
             data.obs_freq_avg = mean(obsavg);

             output = [data output];
             save('spectral_slope_output.mat','output');


        end
    end

end

save('spectral_slope_output.mat','output');