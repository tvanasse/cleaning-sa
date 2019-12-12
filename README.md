## Procssing Pipeline <br />
(1) Extract 5 minutes before awakening <br />
(2) Find that awakening matches to experiment entry <br />
(3) Convert from 256 to "inside" 185 channels <br />
(4) High-pass Filter (1 Hz) <br />
(5) Clean Line Noise <br />
(6) Mark if awakening (30 seconds before) was during N2 or N3 sleep. <br />
(7) Merge NREM awakenings <br />
<br />
(8) Manually remove bad channels from merged file <br />
<br />
(9) Run AMICA on EEG data <br /> 
<br />
(10) Identify artifact components with help of ICLabel <br />
<br />
(11) Split awakenings, and clean bad time streches <br />
(12) High-pass filter (50 Hz) <br />
<br />
(13) Re-merge NREM awakenings and run AMICA again <br />
<br />
(14) Remove artifactual components, and split data again 
