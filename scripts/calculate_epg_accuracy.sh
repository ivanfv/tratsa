#/bin/bash

benchs=("epg.txt")
win="256"
exp="2 3 4 5 6 7 8"
man="2 5 7 10 13 16 20 23"
til="1024 65536 262144 2100000"

for i in "${benchs[@]}"; do
  for e in $exp; do
    for m in $man; do
      for t in $til; do
            	echo "## $i $e $m"
            	./accuracy/check_accuracy references/epg/MP_epg_scaled_w16384_t${t}_double.txt references/epg/MP_epg_scaled_w16384_t${t}_float.txt outputs/epg/MP_epg_scaled.txt_16384_${e}_${m}_${t}.txt ./accuracy_results/epg/accuracy_epg_w16384_e${e}_m${m}_t${t}.txt &
      done;
    done;
  done;
done;

