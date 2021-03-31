# Spectral Exponent


this code allows to compute the spectral exponent of the resting EEG, based on the Power Spectral Density (PSD), over a given scaling region.

The spectral exponent describes the decay of the PSD. it is computed as the slope of an OLS line, fit on log-freq vs log-PSD, excluding oscillatory peaks (and their base).

![figure 1 fit steps](https://user-images.githubusercontent.com/6671316/49792598-72717e00-fd33-11e8-993c-0430313dee43.png)

## USAGE EXAMPLE:
````matlab
% you should have in your workspace
% sRate: sampling Rate, e.g. 1450 for Nexstim
% myEEGch: a vector of datapoints.  

% here, we generate 5 minutes of dummy data, with a 1/f^2 decay, alpha and beta oscillations
sRate= 1450;
dataPoints= 1: sRate*5*60 +2*sRate;
myEEGch=  randn(1,dataPoints(end)); 
myEEGch= smooth( myEEGch, sRate*2)'; 
for ff= [ 9:0.05:11  18:0.05:24] 
    myEEGch= myEEGch+ sin( 2*pi*ff/sRate.*dataPoints)/20/(ff^2); 
end 
myEEGch= myEEGch(sRate+1:end-sRate);


% first compute the PSD
 epLen= 2* sRate; epShift= 1*sRate;numFFT=[];
 [myPSD,frex]= pwelch( myEEGch  , epLen, epShift,numFFT, sRate); 
 
 frBand=[1 40];
 frBins= dsearchn( frex, frBand(1) ):  dsearchn( frex, frBand(2));
 XX= frex(frBins);
 YY= myPSD(frBins);
 robRegMeth= 'ols'; % method to perform linear regression. see >> help robustfit
 doPlot= 1; figure;
 thisCol= [0 0 1];
 [intSlo, stat, Pows, Deviants,  stat0, intSlo0] = fitPowerLaw3steps(XX,YY, robRegMeth,  doPlot, thisCol)
 spectralExponent= intSlo(2);
 
 % repeat for every electrode to obtain the average spectral exponent across the scalp

%% fitPowerLaw3steps %%

%%%%%%%%% INPUT 
% XX ---> are the frequency bins of interest (hz)
% YY ---> are the PSD relative to the bins of interest (mV^2 /Hz)
% robRegMeth --->is the method for robust regression fit (default: ols)
%                type help robustfit for more details
% doPlot  --->    0 (default) means no plot ; 
%                 1 is plot XX and YY on log-log scale;  
%                 2 is plot log(XX) and log(YY) on linear scale; 
% thisCol --->    if plotting, is the rgb triplet
%%%%%%%%% OUTPUT
% intSlo(1) intercept of 2nd (final) powerLaw Fit
% intSlo(2) slope of 2nd (final) powerLaw Fit. i.e. the spectral exponent
% stat is a structure with informations on the 2nd fit, (see 2nd argument of robustfit)

% Pows is a structure with
%           pred: 1*length(frBins) doubles, predicted PSD by the power-law fit
%            obs: 1*length(frBins) doubles, observed PSD    = YY
%            res: 1*length(frBins) doubles, residuals (pred-obs) PSD
%           frex: 1*length(frBins) doubles, observed frex   = XX
%       meanPred: 1 double, mean of the predicted PSD
%        meanObs: 1 double, mean of the observed PSD
%        meanRes: 1 double, mean of the residuals of the PSD
%     meanResPos: 1 double, mean of the residuals of the PSD (exclude negative values)
%       cenFrObs: 1 double, central Frequency of the observed PSD
%       cenFrRes: 1 double, central Frequency of the residual PSD

% Deviants is a structure with 
%       res2: 1*length(frBins) doubles, 
%        res: 1*length(frBins) doubles, 
%        rej:  length(frBins)*1 logicals, 1 if bins are rejected, 0 if kept   for fitting the 2nd power-law
%         fr: min and max freqs considered
%       thre: threshold used for bins adjacent to the peaks
%     clusts: structure with information on the contiguous freq. bins 

% stat0 is a structure with informations on the 1st fit, (see 2nd argument of robustfit)

% intSlo(1) intercept of 2nd (final) powerLaw Fit
% intSlo(2) slope of 2nd (final) powerLaw Fit
````

## HOW DOES IT WORK?
first transform PSD (YY) and frequencies (XX) in log-log and upsample them by 4 times

fit a 1st line, find all the residual >0

small peaks (whose residuals are smaller than 1*mad(residuals)) do not count as peaks

find clusters of rejected frequency bins 

consider only clusters of frequency bins where there is a peak

i.e. reject large enough peaks and their base (adjacent residuals >0)

fit a second line, on the remaining residuals (those closer to a power-law)

the slope of the second line is the spectral exponent. 

get statistics relative to the fit
