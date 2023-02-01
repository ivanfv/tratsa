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
#"audio-MPIII-SVD.txt 200 256 1"S
#"seismology-MPIII-SVE_n180000.txt 50 256 0.01"
#"e0103_n180000.txt 500 256 1"

benchs=("power.txt")
win="256"
exp="2 3 4 5 6 7 8"
man="2 5 7 10 13 16 20 23"
til="1024 65536 262144 2100000"

for i in "${benchs[@]}"; do
  for e in $exp; do
    for m in $man; do
              for t in $til; do
      	echo "## $i $e $m"
	./accuracy/check_accuracy references/power/MP_power_scaled_w1536_t${t}_double.txt outputs/power/MP_power_scaled.txt_1536_8_23_${t}.txt outputs/power/MP_power_scaled.txt_1536_${e}_${m}_${t}.txt ./accuracy_results/power/accuracy_power_w1536_e${e}_m${m}_t${t}.txt &
done;
    done;
  done;
done;

