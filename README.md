## TO-DO <br />
(1) Do multi-taper time/freq <br />

## Processing Pipeline <br />
(1) Extract 5 minutes before awakening <br />
(1) Find that awakening matches to experiment entry <br />
(1) Convert from 256 to "inside" 185 channels <br />
(1) High-pass Filter (1 Hz) <br />
(1) Clean Line Noise <br />
(1) Mark if awakening (0-30 seconds before) was during N2 or N3 sleep. Also, 5 minutes of awakening could not contain any REM to be marked as NREM awakening <br />
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
