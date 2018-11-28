#!/bin/sh
for sub in IR74 IR65 IR71 IR68 IR66 IR63 IR61 CP24 IR67 IR69 IR74 IR78 IR75 IR79
do
cat *log | grep  "Clip: movie" | awk '{print $6}' |  cut  -c 10-13 >> ${outpath}/${sub}_events.txt  # grabbing the movie
cat *log | grep  "Clip: autoLog = True" | awk '{print $1}'  >> ${outpath}/${sub}_timing.txt # grabbing the timing info
paste -d ',' ${outpath}/${sub}_timing.txt ${outpath}/${sub}_events.txt >> ${outpath}/${sub}_eventInfo.txt # combining them
done
