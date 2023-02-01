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

benchs=("penguin_sample_TutorialMPweb.txt 2000 64 1"
        "human_activity-MPIII-SVC.txt 120 64 1")
exp="2 3 4 5 6 7 8"
man="10 12 14 16 18 20 22"

for i in "${benchs[@]}"; do
  for e in $exp; do
    for m in $man; do
      echo "## $i $e $m"
      ./scrimp/scrimp_topK $i $e $m
    done;
  done;
done;

