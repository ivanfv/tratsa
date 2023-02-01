#/bin/bash

benchs=("e0103.txt")
win="256"
exp="2 3 4 5 6 7 8"
man="2 5 7 10 13 16 20 23"
til="1024 65536 262144 2100000"

for i in "${benchs[@]}"; do
  for e in $exp; do
    for m in $man; do
      for t in $til; do
      	echo "## $i $e $m"
	   ./accuracy/check_accuracy references/e0103/MP_e0103_w512_t${t}_double.txt references/e0103/MP_e0103_w512_t${t}_float.txt outputs/e0103/MP_e0103.txt_512_${e}_${m}_${t}.txt ./accuracy_results/e0103/accuracy_e0103_w512_e${e}_m${m}_t${t}.txt &
      done;
    done;
  done;
done;

