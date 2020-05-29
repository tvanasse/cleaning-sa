## TO-DO <br />
(1) Do multi-taper time/freq <br />

## Processing Pipeline <br />
(1) Extract 5 minutes before awakening <br />
(1) Find that awakening matches to experiment entry <br />
(1) Convert from 256 to "inside" 185 channels <br />
(1) High-pass Filter (1 Hz) <br />
(1) Clean Line Noise <br />
(1) Mark if awakening (0-60 seconds before) was during N2 or N3 sleep. Also, 5 minutes of awakening could not contain any REM to be marked as NREM awakening <br />
(1) Merge NREM awakenings <br />
<br />
(2) Manually remove bad channels from merged file <br />
<br />
(3) Run AMICA on EEG data <br /> 
<br />
(4) Identify bad stretches according to component activations, and save rejected data <br />
<br />
(5) Run AMICA on EEG w/ bad streches removed <br />
<br />
(6) Remove artifactual components with help from ICLabel <br />
<br />
(7) Low-pass filter (50 Hz), then re-split awakenings, and perform qa plots.

## Intermediary Data Naming <br />
*_eeg.set - 5 min. raw eeg data before awakening din <br />
*_eeg_hp_trim.set - data high-pass filetered and trimmed (1 sec. beginning and end) <br />
*_eeg_hp_trim_merged.set - merged nrem data <br />
*_eeg_hp_trim_merged_nobadch.set - merged nrem data with bad channels removed <br />
*_eeg_hp_trim_meged_bobadch_ica.set - post ica data on merged nrem data <br />
*-posICA1-cleaned.set - awakening data where bad time-stretches are removed based upon component activations <br />
*-badsection.mat - matrix data of stretches removed per awakening <br />
*nrem_merged_ica2.set - merged nrem awakenings after bad stretch removal and ica <br />
*nrem_merged_ica2_.set - merged nrem awakenings after bad stretch removal and component subtration <br />
*awakening-x-cleaned2_nrem.set cleaned data per awakening <br />
