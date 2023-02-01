#/bin/bash

benchs=("audio.txt")
win="256"
exp="2 3 4 5 6 7 8"
man="2 5 7 10 13 16 20 23"
til="65536 262144 2100000"

for i in "${benchs[@]}"; do
  for e in $exp; do
    for m in $man; do
      for t in $til; do
      	echo "## $i $e $m"
	   ./accuracy/check_accuracy references/audio/MP_audio_quijote_w16384_t${t}_double.txt references/audio/MP_audio_quijote_w16384_t${t}_float.txt outputs/audio/MP_audio.txt_16384_${e}_${m}_${t}.txt ./accuracy_results/audio/accuracy_audio_w16384_e${e}_m${m}_t${t}.txt &
      done;
    done;
  done;
done;

