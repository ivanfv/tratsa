#/bin/bash
#RIC float  exp=8  man=23
#RIC double exp=11 man=52
#usage2: ./scrimp_topK tseries.txt win_size num_threads scale_factor FF_exponent FF_mantissa
#Para usar el HMC y acelerar un poco la cosa:
#numactl --membind 1 ./ejecutable argumentos

#La serie de penguin está sacada del Tutorial de la página https://www.cs.ucr.edu/~eamonn/MatrixProfile.html
# El tamaño de ventana está sacado del paper Matrix Profile II
#./scrimp/scrimp_topK penguin_sample_TutorialMPweb.txt 2000 256 1 5 12

benchs=("power-MPIII-SVF_n180000.txt 1325 256 1"
        "audio-MPIII-SVD.txt 200 256 1"
        "seismology-MPIII-SVE_n180000.txt 50 256 1"
        "e0103_n180000.txt 500 256 1"
        "penguin_sample_TutorialMPweb.txt 800 256 1"
        "human_activity-MPIII-SVC.txt 120 256 1")

for i in "${benchs[@]}"; do
  cfg=`sed 's/txt.*/cfg/g' <<< "$i"`
  #Machaco el archivo y meto exp y man para el cálculo de la distance
  echo "8 7" > ./configs/$cfg
  #Añado al final del archivo el exp y man para el cálculo del dotp
  echo "8 23" >> ./configs/$cfg
  echo "8 7" >> ./configs/$cfg
  echo "8 7" >> ./configs/$cfg
  echo "## $i ## he: 8 hm: 23 | le: 8 lm: 7 ##"
  ./scrimp/scrimp_topK $i

  #Siguiente
  echo "5 10" > ./configs/$cfg
  echo "8 23" >> ./configs/$cfg
  echo "5 10" >> ./configs/$cfg
  echo "5 10" >> ./configs/$cfg
  echo "## $i ## he: 8 hm: 23 | le: 5 lm: 10 ##"
  ./scrimp/scrimp_topK $i

  #Siguiente
  echo "5 2" > ./configs/$cfg
  echo "8 23" >> ./configs/$cfg
  echo "5 2" >> ./configs/$cfg
  echo "5 2" >> ./configs/$cfg
  echo "## $i ## he: 8 hm: 23 | le: 5 lm: 2 ##"
  ./scrimp/scrimp_topK $i

  #Siguiente
  echo "5 10" > ./configs/$cfg
  echo "8 7" >> ./configs/$cfg
  echo "5 10" >> ./configs/$cfg
  echo "5 10" >> ./configs/$cfg
  echo "## $i ## he: 8 hm: 7 | le: 5 lm: 10 ##"
  ./scrimp/scrimp_topK $i

  #Siguiente
  echo "5 2" > ./configs/$cfg
  echo "8 7" >> ./configs/$cfg
  echo "5 2" >> ./configs/$cfg
  echo "5 2" >> ./configs/$cfg
  echo "## $i ## he: 8 hm: 7 | le: 5 lm: 2 ##"
  ./scrimp/scrimp_topK $i
done;

