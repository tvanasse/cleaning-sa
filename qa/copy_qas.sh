cd /Volumes/data/NCCAM3/SA/wDreamReport/aligned/extraction_TJV/
for f in */*/*/qa_figs.tif; do
    cp -v "$f" /Volumes/data/NCCAM3/SA/wDreamReport/aligned/extraction_TJV/scripts/qa/"${f//\//_}"
done

for f in */*/*/log.txt; do
    cp -v "$f" /Volumes/data/NCCAM3/SA/wDreamReport/aligned/extraction_TJV/scripts/qa/"${f//\//_}"
done
