#/bin/bash
#RIC float  exp=8  man=23
#RIC double exp=11 man=52
#usage2: ./scrimp_topK tseries.txt win_size num_threads scale_factor FF_exponent FF_mantissa
#Para usar el HMC y acelerar un poco la cosa:
#numactl --membind 1 ./ejecutable argumentos

#La serie de penguin está sacada del Tutorial de la página https://www.cs.ucr.edu/~eamonn/MatrixProfile.html
# El tamaño de ventana está sacado del paper Matrix Profile II
#./scrimp/scrimp_topK penguin_sample_TutorialMPweb.txt 2000 256 1 5 12

#"power-MPIII-SVF_n180000.txt 1325 256 0.1"
#"audio-MPIII-SVD.txt 200 256 1"
#"seismology-MPIII-SVE_n180000.txt 50 256 0.01"
#"e0103_n180000.txt 500 256 1"

benchs=("seismology_scaled.txt")
win="64"
exp="2 3 4 5 6 7 8"
#exp="6"
#man="20"
man="2 5 7 10 13 16 20 23"
til="65536 262144 2100000"

for i in "${benchs[@]}"; do
  for e in $exp; do
    for m in $man; do
    	for t in $til; do
      		echo "## $i $e $m"
		./exp_${e}/man_${m}/scamp exp_${e}/man_${m}/SCAMP.xclbin timeseries/${i} ${win} ${t} outputs/seismology/MP_${i}_${win}_${e}_${m}_${t}.txt
      	done;
    done;
  done;
done;

