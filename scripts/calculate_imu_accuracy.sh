#/bin/bash


benchs=("imu.txt")
win="256"
exp="2 3 4 5 6 7 8"
man="2 5 7 10 13 16 20 23"
til="1024 65536 262144 2100000"

for i in "${benchs[@]}"; do
  for e in $exp; do
    for m in $man; do
              for t in $til; do
                    	echo "## $i $e $m"
                    	./accuracy/check_accuracy references/imu/MP_imu_w256_t${t}_double.txt references/imu/MP_imu_w256_t${t}_float.txt outputs/imu/MP_imu.txt_256_${e}_${m}_${t}.txt ./accuracy_results/imu/accuracy_imu_w256_e${e}_m${m}_t${t}.txt &
    done;
  done;
done;
done;
