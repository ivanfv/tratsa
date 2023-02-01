#/bin/bash
#RIC float  exp=8  man=23
#RIC double exp=11 man=52
#usage1: ./scamp tseries.txt win_size num_threads scale_factor error_diffusion
#usage2: ./scamp tseries.txt win_size num_threads scale_factor error_diffusion FF_exponent FF_mantissa
#Para usar el HMC y acelerar un poco la cosa:
#numactl --membind 1 ./ejecutable argumentos

#La serie de penguin está sacada del Tutorial de la página https://www.cs.ucr.edu/~eamonn/MatrixProfile.html
# El tamaño de ventana está sacado del paper Matrix Profile II
#./scrimp/scrimp_topK penguin_sample_TutorialMPweb.txt 2000 256 1 5 12

#benchs=("power-MPIII-SVF_n180000.txt 1325 256 1 0"
#        "audio-MPIII-SVD.txt 200 256 1 0"
#        "seismology-MPIII-SVE_n180000.txt 50 256 1 0"
#        "e0103_n180000.txt 500 256 1 0"
#        "penguin_sample_TutorialMPweb.txt 800 256 1 0"
#        "human_activity-MPIII-SVC.txt 120 256 1 0")
        
benchs=("human_activity-MPIII-SVC.txt 120 256 1 0")

exp="2 3 4 5  6  7  8"
man="2 5 7 10 13 16 20 23"

for i in "${benchs[@]}"; do
  for e in $exp; do
    for m in $man; do
      echo "## $i $e $m"
      ./scamp $i $e $m
    done;
  done;
done;

