/* #############################################################################

                 ███████╗ ██████╗██████╗ ██╗███╗   ███╗██████
                 ██╔════╝██╔════╝██╔══██╗██║████╗ ████║██╔══██╗
                 ███████╗██║     ██████╔╝██║██╔████╔██║██████╔╝
                 ╚════██║██║     ██╔══██╗██║██║╚██╔╝██║██╔═══╝ 
                 ███████║╚██████╗██║  ██║██║██║ ╚═╝ ██║██║ 
                 ╚══════╝ ╚═════╝╚═╝  ╚═╝╚═╝╚═╝     ╚═╝╚═╝     
                 +-+-+-+-+-+-+-+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-+
                 | F l e x F l o a t        B e n c h m a r k |
                 +-+-+-+-+-+-+-+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-+

SCRIMP FlexFloat benchmark provides a precision explotation tool for time series
analysis using Matrix Profile. This program allows changing the number of bits
used by floating point operations to meassure how affects to the accuracy of 
result. This program is not inteded to be used for performance evaluation
purposes as long as it is based on a floating point emulation library.

The code is based on the SCRIMP C++ implementation developed by Yan Zhu, 
Chin-Chia Michael Yeh, Zachary Zimmerman, Kaveh Kamgar and Eamonn Keogh that 
can be found at https://sites.google.com/site/scrimpplusplus/. It was modified
by Ivan Fernandez in August 2019 to include the FlexFloat library.

Details of the SCRIMP++ algorithm can be found at:
    - "SCRIMP++: Motif Discovery at Interactive Speeds", ICDM 2018.
    - Matrix Profile website: https://www.cs.ucr.edu/~eamonn/MatrixProfile.html

Details of Flexfloat library can be found at:
    - "FlexFloat: A Software Library for Transprecision Computing",
             IEEE Transactions on Computer-Aided Design of Int. Cir. and Systems
    - FlexFloat repository: https://github.com/oprecomp/flexfloat    

 ********************************************************************************
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
 ********************************************************************************
Usage: 
>> scrimpplusplus InputFile WindowSize nThreads Scale
        - InputFile: Name of the time series file
        - WindowSize: Subsequence length m
        - nThreads: Number of threads to be spawn
        - Scale: Scale factor for the time series data

Example input:
>> scrimp_ff taxi.txt 100 8 1

Example output:
The code will generate one output file (e.g. result_taxi.csv)
        - First column is the index of the data.
        - Second column is the matrix profile flexfloat distance value.
        - Third column is the matrix profile flexfloat index.
        - Fourth column is the matrix profile distance value.
        - Fifth column is the matrix profile index.
        - Sixth column is the error
############################################################################# */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include <unistd.h> //RIC para getpid(), usado para obtener un nombre de fichero único
#include <assert.h> //Meto el assert
#include <signal.h> //Para atrapar señales de fenv.h
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <omp.h>
#include "../flexfloat/include/flexfloat.h"

//RIC macro de valor absoluto para antes del sqrt en algunos valores negativos
//que da flexfloat al calcular las distancias
#define ABS_FF(val,prec) {                                                 \
  flexfloat_t zero_tmp;                                                    \
  ff_init_double(&zero_tmp, 0.0, prec);                                    \
  if (ff_lt(val,&zero_tmp,tid)) ff_inverse(val,val,tid);                   \
  }
//En lugar de hacer el valor absoluto podemos probar poniéndo la distancia a 0
#define ZERO_FF(val,prec) {                                                \
  flexfloat_t zero_tmp;                                                    \
  ff_init_double(&zero_tmp, 0.0, prec);                                    \
  if (ff_lt(val,&zero_tmp,tid)) ff_init_double(val, 0.0, prec);            \
  }

#define PATH_TSERIES "./timeseries/"
#define PATH_CFG "./configs/"
#define PATH_RESULT "./results/result_"
#define EXCLUSION_FACTOR 4
#define RANGE 10 //RIC rango para calcular Top-K accuracy
#define K 100    //RIC top 100

//RIC macros para manejar las excepciones de punto flotante
int *exceptions;
#define INITEXC(numTh)  {                                 \
    exceptions = (int *) malloc(sizeof (int)*numTh * 5);  \
    assert(exceptions);                                   \
  }
#define RESETEXC(numTh) {                                 \
    for(int i=0; i<numTh*5; i++) exceptions[i] = 0;       \
    feclearexcept(FE_ALL_EXCEPT);                         \
  }
#define SAVEEXC(th)   {                       \
    if(fetestexcept(FE_DIVBYZERO)) exceptions[th*5+0]++;  \
    if(fetestexcept(FE_INEXACT)) exceptions[th*5+1]++;    \
    if(fetestexcept(FE_INVALID)) exceptions[th*5+2]++;    \
    if(fetestexcept(FE_OVERFLOW)) exceptions[th*5+3]++;   \
    if(fetestexcept(FE_UNDERFLOW)) exceptions[th*5+4]++;  \
  }
#define PRINTEXC(numTh)  {                         \
    for(int i=5; i<numTh*5; i+=5) {                \
      exceptions[0]+= exceptions[i];                \
      exceptions[1]+= exceptions[i+1];              \
      exceptions[2]+= exceptions[i+2];              \
      exceptions[3]+= exceptions[i+3];              \
      exceptions[4]+= exceptions[i+4];              \
    }                                              \
    printf("[INFO] DivByZero(%d)Inexact(%d)Invalid(%d)Overflow(%d)Underflow(%d)\n",exceptions[0],exceptions[1],exceptions[2],exceptions[3],exceptions[4]); \
  }


//RIC para hacer la ordenación del perfil necesito una estructura que contenga
//el valor del perfil y los dos índices (uno para la ventana i y otro para la j)
//ya que al ordenar el mprofile se pierde la información de la i

typedef struct ProfItem_f {
  float val;
  int i;
  int j;
} profItem_f;

typedef struct ProfItem_d {
  double val;
  int i;
  int j;
} profItem_d;

typedef struct ProfItem_ff {
  flexfloat_t val;
  int i;
  int j;
} profItem_ff;

unsigned dist_exp, dist_man;
unsigned prof_exp, prof_man;
unsigned dotp_exp, dotp_man;
unsigned stats_exp, stats_man;

//RIC función de comparación para qsort

int compar_f(const void * a, const void * b) {
  float aval, bval;
  aval = ((profItem_f *) a)->val;
  bval = ((profItem_f *) b)->val;
  if (aval < bval) return -1;
  if (aval == bval) return 0;
  if (aval > bval) return 1;
}

int compar_d(const void * a, const void * b) {
  double aval, bval;
  aval = ((profItem_d *) a)->val;
  bval = ((profItem_d *) b)->val;
  if (aval < bval) return -1;
  if (aval == bval) return 0;
  if (aval > bval) return 1;
}

int compar_ff(const void * a, const void * b) {
  flexfloat_t aval, bval;
  aval = ((profItem_ff *) a)->val;
  bval = ((profItem_ff *) b)->val;
  if (ff_lt(&aval, &bval, 0)) return -1;
  if (ff_eq(&aval, &bval, 0)) return 0;
  if (ff_gt(&aval, &bval, 0)) return 1;
}

void scrimp_f(float * tSeries, float * AMean, float * ASigma, int ProfileLength,
        int windowSize, int* idx, profItem_f * profile_f, int exclusionZone, int numThreads);
void scrimp_d(double * tSeries, double * AMean, double * ASigma, int ProfileLength,
        int windowSize, int* idx, profItem_d * profile_d, int exclusionZone, int numThreads);
void scrimp_ff(flexfloat_t * tSeries, flexfloat_t* AMean, flexfloat_t* ASigma,
        int ProfileLength, int windowSize, int* idx, profItem_ff * profile_ff,
        int exclusionZone, int numThreads);

void print_header();
void read_config(char * file_name);
void getProfileAccuracy(char *result_path, int ProfileLength, profItem_f *profile_f,
        profItem_d *profile_d, profItem_ff *profile_ff);
void ffStatsToFile(char *result_path, int numThreads);

int main(int argc, char* argv[]) {
  FILE * fp;
  char path_tSeries[1000], path_result[1000];
  int windowSize, ProfileLength, exclusionZone, timeSeriesLength, numThreads, * idx;
  float * tSeries_f, * AMean_f, * ASigma_f;
  profItem_f *profile_f;
  double tSeriesMin_d, tSeriesMax_d, scaleFactor_d, * tSeries_d, * AMean_d, * ASigma_d;
  profItem_d *profile_d;
  flexfloat_t * tSeries_ff, * AMean_ff, * ASigma_ff;
  profItem_ff * profile_ff;
  struct timeval tIni, tEnd; //RIC 
  double mp_time_f, mp_time_d, mp_time_ff; //RIC para guardar el tiempo del matrix profile

  print_header();
  /* Getting program arguments ---------------------------------------- */
  if (argc != 5 && argc != 7) {
    printf("[ERROR] usage1: ./scrimp_topK tseries.txt win_size num_threads scale_factor\n");
    printf("[ERROR] usage2: ./scrimp_topK tseries.txt win_size num_threads scale_factor FF_exponent FF_mantissa\n");
    return 0;
  }
  windowSize = atoi(argv[2]);
  numThreads = atoi(argv[3]);
  scaleFactor_d = atof(argv[4]);
  omp_set_num_threads(numThreads);
  if (argc == 5)
    read_config(argv[1]);
  else {//RIC argc == 7 pongo a todos el mismo exponente y mantisa
    unsigned exp_tmp = atoi(argv[5]);
    unsigned man_tmp = atoi(argv[6]);
    dist_exp = exp_tmp;
    dist_man = man_tmp;
    dotp_exp = exp_tmp;
    dotp_man = man_tmp;
    stats_exp = exp_tmp;
    stats_man = man_tmp;
    prof_exp = exp_tmp;
    prof_man = man_tmp;
  }
  /* ------------------------------------------------------------------ */
  INITEXC(numThreads);
  ff_start_stats(numThreads);

  /* Time series loading and memory allocating ------------------------ */
  sprintf(path_tSeries, "%s%s", PATH_TSERIES, argv[1]);
  printf("##############################################\n");
  printf("[INFO] Loading %s ...\n", argv[1]);
  timeSeriesLength = 0;
  fp = fopen(path_tSeries, "r");
  assert(fp);
  for (char c = getc(fp); !feof(fp); c = getc(fp))
    if (c == '\n') timeSeriesLength = timeSeriesLength + 1;
  fclose(fp);

  ProfileLength = timeSeriesLength - windowSize + 1;
  exclusionZone = windowSize / EXCLUSION_FACTOR;

  //RIC reservo memoria para el vector de diagonales
  idx = (int *) malloc(sizeof (int) * ProfileLength);
  assert(idx);
  //RIC reservo memoria para float
  tSeries_f = (float *) malloc(sizeof (float) * timeSeriesLength);
  assert(tSeries_f);
  AMean_f = (float *) malloc(sizeof (float) * ProfileLength);
  assert(AMean_f);
  ASigma_f = (float *) malloc(sizeof (float) * ProfileLength);
  assert(ASigma_f);
  profile_f = (profItem_f *) malloc(sizeof (profItem_f) * ProfileLength);
  assert(profile_f);
  //RIC reservo memoria para double
  tSeries_d = (double *) malloc(sizeof (double) * timeSeriesLength);
  assert(tSeries_d);
  AMean_d = (double *) malloc(sizeof (double) * ProfileLength);
  assert(AMean_d);
  ASigma_d = (double *) malloc(sizeof (double) * ProfileLength);
  assert(ASigma_d);
  profile_d = (profItem_d *) malloc(sizeof (profItem_d) * ProfileLength);
  assert(profile_d);
  //RIC reservo memoria para flexfloat
  tSeries_ff = (flexfloat_t *) malloc(sizeof (flexfloat_t) * timeSeriesLength);
  assert(tSeries_ff);
  AMean_ff = (flexfloat_t *) malloc(sizeof (flexfloat_t) * ProfileLength);
  assert(AMean_ff);
  ASigma_ff = (flexfloat_t *) malloc(sizeof (flexfloat_t) * ProfileLength);
  assert(ASigma_ff);
  profile_ff = (profItem_ff *) malloc(sizeof (profItem_ff) * ProfileLength);
  assert(profile_ff);

  //RIC Los mínimos y máximos de la serie se pueden iniciar al primer valor de la misma
  //RIC cambio tSeriesMax = 0;
  tSeriesMin_d = INFINITY;
  tSeriesMax_d = -INFINITY; //RIC puede ser la serie entera de valores negativos?
  //RIC leo la serie de fichero
  fp = fopen(path_tSeries, "r");
  assert(fp);
  for (int i = 0; i < timeSeriesLength; i++) {
    int err = fscanf(fp, "%lf", &tSeries_d[i]);
    if (err < 0) return -1;
    tSeries_d[i] *= scaleFactor_d;
    if (tSeries_d[i] < tSeriesMin_d) tSeriesMin_d = tSeries_d[i];
    if (tSeries_d[i] > tSeriesMax_d) tSeriesMax_d = tSeries_d[i];

    tSeries_f[i] = (float) tSeries_d[i];
    ff_init_double(&tSeries_ff[i], tSeries_d[i], (flexfloat_desc_t) {dotp_exp, dotp_man});
  }
  fclose(fp);
  printf("[INFO] DONE\n");
  /* ------------------------------------------------------------------ */

  printf("[INFO] Program parameters:\n");
  printf("----------------------------------------------\n");
  printf("  Time series length: %d\n", timeSeriesLength);
  printf("  Window size:        %d\n", windowSize);
  printf("  Time series min:    %f\n", tSeriesMin_d);
  printf("  Time series max:    %f\n", tSeriesMax_d);
  printf("  Dot product max:    %f\n", pow(tSeriesMax_d, 2) * windowSize);
  printf("  Number of threads:  %d\n", numThreads);
  printf("  Scale factor:       %.4f\n", scaleFactor_d);
  printf("  FF dist - exp, man: %d, %d\n", dist_exp, dist_man);
  printf("  FF dotp - exp, man: %d, %d\n", dotp_exp, dotp_man);
  printf("  FF stat - exp, man: %d, %d\n", stats_exp, stats_man);
  printf("  FF prof - exp, man: %d, %d\n", prof_exp, prof_man);
  printf("----------------------------------------------\n");

  /* Preprocessing statistics ----------------------------------------- */
  printf("[INFO] Preprocessing statistics ...\n");
  double * ACumSum = (double *) malloc(sizeof (double) * timeSeriesLength);
  assert(ACumSum);
  double * ASum = (double *) malloc(sizeof (double) * timeSeriesLength);
  assert(ASum);
  double * ASumSq = (double *) malloc(sizeof (double) * timeSeriesLength);
  assert(ASumSq);
  double * ASigmaSq = (double *) malloc(sizeof (double) * timeSeriesLength);
  assert(ASigmaSq);
  double * ASqCumSum = (double *) malloc(sizeof (double) * timeSeriesLength);
  assert(ASqCumSum);

  //RIC para las excepciones de punto flotante. Las limpio antes de los cáculos
  RESETEXC(numThreads);
  //RIC suma acumulada. Cada elemento de ACumSum contiene la suma de los anteriores
  ACumSum[0] = tSeries_d[0];
  //RIC lo mismo pero con la suma de cuadrados para calcular la desviación típica
  ASqCumSum[0] = tSeries_d[0] * tSeries_d[0];
  for (int i = 1; i < timeSeriesLength; i++) {
    ACumSum[i] = tSeries_d[i] + ACumSum[i - 1];
    ASqCumSum[i] = tSeries_d[i] * tSeries_d[i] + ASqCumSum[i - 1];
  }
  //RIC por cada ventana se resta al siguiente del último el de antes del primero.
  //    Así se obtiene la suma de todos los elementos de la ventana
  ASum[0] = ACumSum[windowSize - 1];
  //RIC lo mismo para la suma de cuadrados
  ASumSq[0] = ASqCumSum[windowSize - 1];
  for (int i = 0; i < timeSeriesLength - windowSize; i++) {
    ASum[i + 1] = ACumSum[windowSize + i] - ACumSum[i];
    ASumSq[i + 1] = ASqCumSum[windowSize + i] - ASqCumSum[i];
  }
  //RIC dividir por el tamaño de ventana para obtener la media
  for (int i = 0; i < ProfileLength; i++) {
    AMean_d[i] = ASum[i] / windowSize;
    //RIC la desviación típica se puede calcular así:
    ASigmaSq[i] = ASumSq[i] / windowSize - AMean_d[i] * AMean_d[i];
    ASigma_d[i] = sqrt(ASigmaSq[i]);

    //RIC aprovecho para inicializar profiles y otros valores para otras precisiones
    //RIC no hace falta iniciar el profile pq se machaca en scrimp
    //RIC además parece que iniciar con INFINITY una variable flexfloat hace settear el bit de overflow
    //profile_f[i].val = INFINITY;
    //profile_d[i].val = INFINITY;
    //ff_init_double(&(profile_ff[i].val), profile_d[i].val, (flexfloat_desc_t) {prof_exp, prof_man});
    AMean_f[i] = (float) AMean_d[i];
    ASigma_f[i] = (float) ASigma_d[i];

    ff_init_double(&AMean_ff[i], AMean_d[i], (flexfloat_desc_t) {stats_exp, stats_man});
    ff_init_double(&ASigma_ff[i], ASigma_d[i], (flexfloat_desc_t) {stats_exp, stats_man});
  }
  free(ACumSum);
  free(ASum);
  free(ASqCumSum);
  free(ASumSq);
  free(ASigmaSq);
  //RIC después de hacer los cálculos miro si ha habido excepciones
  SAVEEXC(0);
  PRINTEXC(numThreads);
  printf("[INFO] DONE\n");
  /* ------------------------------------------------------------------ */
  //RIC le quito el .txt a argv[1]
  argv[1][strlen(argv[1]) - 4] = 0;
  //RIC creo el nombre del fichero con sprintf y meto el pid para crear un nombre único  
  sprintf(path_result, "%s%s_w%d_s%.2f_t%d_dotp%d-%d_dist%d-%d_stat%d-%d_prof%d-%d_%d.csv", PATH_RESULT, argv[1], windowSize,
          scaleFactor_d, numThreads, dotp_exp, dotp_man, dist_exp, dist_man, stats_exp, stats_man, prof_exp, prof_man, getpid());
  fp = fopen(path_result, "w");
  assert(fp);
  fprintf(fp, "#Stats Exceptions: DivByZero,Inexact,Invalid,Overflow,Underflow\n");
  fprintf(fp, "%d,%d,%d,%d,%d\n", exceptions[0], exceptions[1], exceptions[2], exceptions[3], exceptions[4]);

  /* Initializing diagonals and flexfloat statistics------------------- */
  for (int i = exclusionZone + 1; i < ProfileLength; i++)
    idx[i - (exclusionZone + 1)] = i;
  /* ------------------------------------------------------------------ */
  // Running SCRIMP double precision-----------------------------------
  printf("[INFO] Running SCRIMP double ...\n");
  RESETEXC(numThreads);
  gettimeofday(&tIni, NULL);
  scrimp_d(tSeries_d, AMean_d, ASigma_d, ProfileLength, windowSize, idx,
          profile_d, exclusionZone, numThreads);
  gettimeofday(&tEnd, NULL);
  PRINTEXC(numThreads);
  fprintf(fp, "#Scrimp Double Exceptions: DivByZero,Inexact,Invalid,Overflow,Underflow\n");
  fprintf(fp, "%d,%d,%d,%d,%d\n", exceptions[0], exceptions[1], exceptions[2], exceptions[3], exceptions[4]);
  mp_time_d = (double) (tEnd.tv_sec - tIni.tv_sec) + (double) (tEnd.tv_usec - tIni.tv_usec) / 1000000;
  printf("[INFO] DONE (elapsed time %.3f seconds)\n", mp_time_d);
  /* ------------------------------------------------------------------ */
  // Running SCRIMP single precision ----------------------------------
  printf("[INFO] Running SCRIMP float ...\n");
  RESETEXC(numThreads);
  gettimeofday(&tIni, NULL);
  scrimp_f(tSeries_f, AMean_f, ASigma_f, ProfileLength, windowSize, idx,
          profile_f, exclusionZone, numThreads);
  gettimeofday(&tEnd, NULL);
  PRINTEXC(numThreads);
  fprintf(fp, "#Scrimp Float Exceptions: DivByZero,Inexact,Invalid,Overflow,Underflow\n");
  fprintf(fp, "%d,%d,%d,%d,%d\n", exceptions[0], exceptions[1], exceptions[2], exceptions[3], exceptions[4]);
  mp_time_f = (double) (tEnd.tv_sec - tIni.tv_sec) + (double) (tEnd.tv_usec - tIni.tv_usec) / 1000000;
  printf("[INFO] DONE (elapsed time %.3f seconds)\n", mp_time_f);
  /* ------------------------------------------------------------------ */
  // Running SCRIMP FlexFloat -----------------------------------------
  printf("[INFO] Running SCRIMP FlexFloat ...\n");
  RESETEXC(numThreads);
  gettimeofday(&tIni, NULL);
  scrimp_ff(tSeries_ff, AMean_ff, ASigma_ff, ProfileLength, windowSize, idx,
          profile_ff, exclusionZone, numThreads);
  gettimeofday(&tEnd, NULL);
  PRINTEXC(numThreads);
  fprintf(fp, "#Scrimp FlexFloat Exceptions: DivByZero,Inexact,Invalid,Overflow,Underflow\n");
  fprintf(fp, "%d,%d,%d,%d,%d\n", exceptions[0], exceptions[1], exceptions[2], exceptions[3], exceptions[4]);
  mp_time_ff = (double) (tEnd.tv_sec - tIni.tv_sec) + (double) (tEnd.tv_usec - tIni.tv_usec) / 1000000;
  printf("[INFO] DONE (elapsed time %.3f seconds)\n", mp_time_ff);

  /* ------------------------------------------------------------------ */
  // Saving results to file ------------------------------------------- 
  printf("[INFO] Saving %s ...\n", path_result);

  //RIC imprimo el tiempo en primer lugar
  fprintf(fp, "#Time (s): double,float,flexfloat\n");
  fprintf(fp, "%.9f,%.9f,%9f\n", mp_time_d, mp_time_f, mp_time_ff);
  fprintf(fp, "#Profile Length\n");
  fprintf(fp, "%d\n", ProfileLength);
  fprintf(fp, "#i,tseries,doubleMP,j,floatMP,j,errorf(%%),flexfloatMP,j,errorff(%%)\n");
  double errorf, errorff;
  for (int i = 0; i < ProfileLength; i++) {
    errorf = (fabs((sqrt(profile_d[i].val) - sqrt(profile_f[i].val))) / sqrt(profile_d[i].val)) * 100;
    errorff = (fabs((sqrt(profile_d[i].val) - sqrt(ff_get_double(&(profile_ff[i].val))))) / sqrt(profile_d[i].val)) * 100;
    fprintf(fp, "%d,%.9f,%.9f,%d,%.9f,%d,%.9f,%.9f,%d,%.9f\n", i, tSeries_d[i], sqrt(profile_d[i].val),
            profile_d[i].i, (float) sqrt(profile_f[i].val), profile_f[i].i, errorf,
            sqrt(ff_get_double(&(profile_ff[i].val))), profile_ff[i].i, errorff);
  }
  fclose(fp);
  /* ------------------------------------------------------------------ */
  printf("[INFO] FlexFloat stats:\n");
  printf("----------------------------------------------\n");
  ff_print_stats();
  ffStatsToFile(path_result, numThreads);
  printf("----------------------------------------------\n");
  // ------------------------------------------------------------------
  // RIC Qsorting Mprofiles and extracting top K motifs and discords accuracy
  getProfileAccuracy(path_result, ProfileLength, profile_f, profile_d, profile_ff);

  printf("[INFO] DONE\n");
  /* ------------------------------------------------------------------ */
  free(tSeries_f);
  free(AMean_f);
  free(ASigma_f);
  free(profile_f);
  free(tSeries_d);
  free(AMean_d);
  free(ASigma_d);
  free(profile_d);
  free(tSeries_ff);
  free(AMean_ff);
  free(ASigma_ff);
  free(profile_ff);
  free(idx);
  free(exceptions);
  printf("##############################################\n");

  return 0;
}

void ffStatsToFile(char *result_path, int numThreads) {
  /*--------------------------------------------------------------------------*/
  // RIC imprimo al final del fichero
  FILE *fp = fopen(result_path, "a");
  assert(fp);
  OpStats statsAcc;
  statsAcc.minus = 0;
  statsAcc.add = 0;
  statsAcc.sub = 0;
  statsAcc.mul = 0;
  statsAcc.div = 0;
  statsAcc.fma = 0;
  statsAcc.cmp = 0;
  for (int i = 0; i < numThreads; i++) {
    OpStats * stats = getOpStats((flexfloat_desc_t) {dist_exp, dist_man}, i);
    statsAcc.minus += stats->minus;
    statsAcc.add += stats->add;
    statsAcc.sub += stats->sub;
    statsAcc.mul += stats->mul;
    statsAcc.div += stats->div;
    statsAcc.fma += stats->fma;
    statsAcc.cmp += stats->cmp;
  }
  fprintf(fp, "#Flexfloat<%d,%d> dist stats: INV,ADD,SUB,MUL,DIV,FMA,CMP\n", dist_exp, dist_man);
  fprintf(fp, "%lu,%lu,%lu,%lu,%lu,%lu,%lu\n", statsAcc.minus, statsAcc.add, statsAcc.sub, statsAcc.mul,
          statsAcc.div, statsAcc.fma, statsAcc.cmp);
  
  statsAcc.minus = 0;
  statsAcc.add = 0;
  statsAcc.sub = 0;
  statsAcc.mul = 0;
  statsAcc.div = 0;
  statsAcc.fma = 0;
  statsAcc.cmp = 0;
  for (int i = 0; i < numThreads; i++) {
    OpStats * stats = getOpStats((flexfloat_desc_t) {dotp_exp, dotp_man}, i);
    statsAcc.minus += stats->minus;
    statsAcc.add += stats->add;
    statsAcc.sub += stats->sub;
    statsAcc.mul += stats->mul;
    statsAcc.div += stats->div;
    statsAcc.fma += stats->fma;
    statsAcc.cmp += stats->cmp;
  }
  fprintf(fp, "#Flexfloat<%d,%d> dotp stats: INV,ADD,SUB,MUL,DIV,FMA,CMP\n", dotp_exp, dotp_man);
  fprintf(fp, "%lu,%lu,%lu,%lu,%lu,%lu,%lu\n", statsAcc.minus, statsAcc.add, statsAcc.sub, statsAcc.mul,
          statsAcc.div, statsAcc.fma, statsAcc.cmp);
  
  statsAcc.minus = 0;
  statsAcc.add = 0;
  statsAcc.sub = 0;
  statsAcc.mul = 0;
  statsAcc.div = 0;
  statsAcc.fma = 0;
  statsAcc.cmp = 0;
  for (int i = 0; i < numThreads; i++) {
    OpStats * stats = getOpStats((flexfloat_desc_t) {stats_exp, stats_man}, i);
    statsAcc.minus += stats->minus;
    statsAcc.add += stats->add;
    statsAcc.sub += stats->sub;
    statsAcc.mul += stats->mul;
    statsAcc.div += stats->div;
    statsAcc.fma += stats->fma;
    statsAcc.cmp += stats->cmp;
  }
  fprintf(fp, "#Flexfloat<%d,%d> stats stats: INV,ADD,SUB,MUL,DIV,FMA,CMP\n", stats_exp, stats_man);
  fprintf(fp, "%lu,%lu,%lu,%lu,%lu,%lu,%lu\n", statsAcc.minus, statsAcc.add, statsAcc.sub, statsAcc.mul,
          statsAcc.div, statsAcc.fma, statsAcc.cmp);

  statsAcc.minus = 0;
  statsAcc.add = 0;
  statsAcc.sub = 0;
  statsAcc.mul = 0;
  statsAcc.div = 0;
  statsAcc.fma = 0;
  statsAcc.cmp = 0;
  for (int i = 0; i < numThreads; i++) {
    OpStats * stats = getOpStats((flexfloat_desc_t) {prof_exp, prof_man}, i);
    statsAcc.minus += stats->minus;
    statsAcc.add += stats->add;
    statsAcc.sub += stats->sub;
    statsAcc.mul += stats->mul;
    statsAcc.div += stats->div;
    statsAcc.fma += stats->fma;
    statsAcc.cmp += stats->cmp;
  }
  fprintf(fp, "#Flexfloat<%d,%d> prof stats: INV,ADD,SUB,MUL,DIV,FMA,CMP\n", prof_exp, prof_man);
  fprintf(fp, "%lu,%lu,%lu,%lu,%lu,%lu,%lu\n", statsAcc.minus, statsAcc.add, statsAcc.sub, statsAcc.mul,
          statsAcc.div, statsAcc.fma, statsAcc.cmp);
  //CASTS
  CastStats castAcc;
  castAcc.total = 0;
  for (int i = 0; i < numThreads; i++) {
    CastStats * stats = getCastStats((flexfloat_desc_t) {dotp_exp, dotp_man},(flexfloat_desc_t) {dist_exp, dist_man}, i);
    castAcc.total += stats->total;
  }
  fprintf(fp, "#Flexfloat<%d,%d> dotp to Flexfloat<%d,%d> dist casts\n", dotp_exp, dotp_man, dist_exp, dist_man);
  fprintf(fp, "%lu\n", castAcc.total);
  
  castAcc.total = 0;
  for (int i = 0; i < numThreads; i++) {
    CastStats * stats = getCastStats((flexfloat_desc_t) {stats_exp, stats_man}, (flexfloat_desc_t) {dist_exp, dist_man}, i);
    castAcc.total += stats->total;
  }
  fprintf(fp, "#Flexfloat<%d,%d> stats to Flexfloat<%d,%d> dist casts\n", stats_exp, stats_man, dist_exp, dist_man);
  fprintf(fp, "%lu\n", castAcc.total);
  
  castAcc.total = 0;
  for (int i = 0; i < numThreads; i++) {
    CastStats * stats = getCastStats((flexfloat_desc_t) {dist_exp, dist_man}, (flexfloat_desc_t) {prof_exp, prof_man}, i);
    castAcc.total += stats->total;
  }
  fprintf(fp, "#Flexfloat<%d,%d> dist to Flexfloat<%d,%d> prof casts\n", dist_exp, dist_man, prof_exp, prof_man);
  fprintf(fp, "%lu\n", castAcc.total);
  
  fclose(fp);
}

void getProfileAccuracy(char *path_result, int ProfileLength, profItem_f *profile_f,
        profItem_d *profile_d, profItem_ff *profile_ff) {
  struct timeval tIni, tEnd; //RIC 
  double qs_time;

  printf("[INFO] QSorting MProfiles ...\n");
  gettimeofday(&tIni, NULL);
  qsort(profile_d, ProfileLength, sizeof (profItem_d), compar_d);
  qsort(profile_f, ProfileLength, sizeof (profItem_f), compar_f);
  qsort(profile_ff, ProfileLength, sizeof (profItem_ff), compar_ff);
  gettimeofday(&tEnd, NULL);
  qs_time = (double) (tEnd.tv_sec - tIni.tv_sec) + (double) (tEnd.tv_usec - tIni.tv_usec) / 1000000;
  printf("[INFO] DONE (elapsed time %.3f seconds)\n", qs_time);

  //RIC calculo el porcentaje de matches del top-K motifs del transpreciso y el float con el double
  int motifs_f = 0, motifs_f_range = 0, motifs_ff = 0, motifs_ff_range = 0;
  for (int i = 0; i < K; i++) {
    for (int j = 0; j < K; j++) {
      //RIC me da igual el orden (i,j) (j,i)
      if (((profile_f[i].i == profile_d[j].i) && (profile_f[i].j == profile_d[j].j)) ||
              ((profile_f[i].j == profile_d[j].i) && (profile_f[i].i == profile_d[j].j))) {
        motifs_f++;
        break;
      }
    }
    for (int j = 0; j < K; j++) {
      if ((((profile_f[i].i <= profile_d[j].i + RANGE) && (profile_f[i].i >= profile_d[j].i - RANGE)) &&
              ((profile_f[i].j <= profile_d[j].j + RANGE) && (profile_f[i].j >= profile_d[j].j - RANGE))) ||
              (((profile_f[i].j <= profile_d[j].i + RANGE) && (profile_f[i].j >= profile_d[j].i - RANGE)) &&
              ((profile_f[i].i <= profile_d[j].j + RANGE) && (profile_f[i].i >= profile_d[j].j - RANGE)))) {
        motifs_f_range++;
        break;
      }
    }
    for (int j = 0; j < K; j++) {
      //RIC me da igual el orden (i,j) (j,i)
      if (((profile_ff[i].i == profile_d[j].i) && (profile_ff[i].j == profile_d[j].j)) ||
              ((profile_ff[i].j == profile_d[j].i) && (profile_ff[i].i == profile_d[j].j))) {
        motifs_ff++;
        break;
      }
    }
    for (int j = 0; j < K; j++) {
      if ((((profile_ff[i].i <= profile_d[j].i + RANGE) && (profile_ff[i].i >= profile_d[j].i - RANGE)) &&
              ((profile_ff[i].j <= profile_d[j].j + RANGE) && (profile_ff[i].j >= profile_d[j].j - RANGE))) ||
              (((profile_ff[i].j <= profile_d[j].i + RANGE) && (profile_ff[i].j >= profile_d[j].i - RANGE)) &&
              ((profile_ff[i].i <= profile_d[j].j + RANGE) && (profile_ff[i].i >= profile_d[j].j - RANGE)))) {
        motifs_ff_range++;
        break;
      }
    }
  }
  //RIC calculo el porcentaje de matches del top-K del transpreciso con el preciso
  int discords_f = 0, discords_f_range = 0, discords_ff = 0, discords_ff_range = 0;
  for (int i = ProfileLength - 1; i >= (ProfileLength - K); i--) {
    for (int j = ProfileLength - 1; j >= (ProfileLength - K); j--) {
      if (((profile_f[i].i == profile_d[j].i) && (profile_f[i].j == profile_d[j].j)) ||
              ((profile_f[i].j == profile_d[j].i) && (profile_f[i].i == profile_d[j].j))) {
        discords_f++;
        break;
      }
    }
    for (int j = ProfileLength - 1; j >= (ProfileLength - K); j--) {
      if ((((profile_f[i].i <= profile_d[j].i + RANGE) && (profile_f[i].i >= profile_d[j].i - RANGE)) &&
              ((profile_f[i].j <= profile_d[j].j + RANGE) && (profile_f[i].j >= profile_d[j].j - RANGE))) ||
              (((profile_f[i].j <= profile_d[j].i + RANGE) && (profile_f[i].j >= profile_d[j].i - RANGE)) &&
              ((profile_f[i].i <= profile_d[j].j + RANGE) && (profile_f[i].i >= profile_d[j].j - RANGE)))) {
        discords_f_range++;
        break;
      }
    }
    for (int j = ProfileLength - 1; j >= (ProfileLength - K); j--) {
      if (((profile_ff[i].i == profile_d[j].i) && (profile_ff[i].j == profile_d[j].j)) ||
              ((profile_ff[i].j == profile_d[j].i) && (profile_ff[i].i == profile_d[j].j))) {
        discords_ff++;
        break;
      }
    }
    for (int j = ProfileLength - 1; j >= (ProfileLength - K); j--) {
      if ((((profile_ff[i].i <= profile_d[j].i + RANGE) && (profile_ff[i].i >= profile_d[j].i - RANGE)) &&
              ((profile_ff[i].j <= profile_d[j].j + RANGE) && (profile_ff[i].j >= profile_d[j].j - RANGE))) ||
              (((profile_ff[i].j <= profile_d[j].i + RANGE) && (profile_ff[i].j >= profile_d[j].i - RANGE)) &&
              ((profile_ff[i].i <= profile_d[j].j + RANGE) && (profile_ff[i].i >= profile_d[j].j - RANGE)))) {
        discords_ff_range++;
        break;
      }
    }
  }
  printf("[INFO] Top-%d Accuracy: \n", K);
  printf("-----------------------------------------------\n");
  printf("Motifs (float):    %.2f%%\n", (double) (motifs_f * 100.0f) / K);
  printf("Motifs (float)(\u00B1%d):    %.2f%%\n", RANGE, (double) (motifs_f_range * 100.0f) / K);
  printf("Discords (float): %.2f%%\n", (double) (discords_f * 100.0f) / K);
  printf("Discords (float)(\u00B1%d): %.2f%%\n\n", RANGE, (double) (discords_f_range * 100.0f) / K);

  printf("Motifs (FlexFloat):    %.2f%%\n", (double) (motifs_ff * 100.0f) / K);
  printf("Motifs (FlexFloat)(\u00B1%d):    %.2f%%\n", RANGE, (double) (motifs_ff_range * 100.0f) / K);
  printf("Discords (FlexFloat): %.2f%%\n", (double) (discords_ff * 100.0f) / K);
  printf("Discords (FlexFloat)(\u00B1%d): %.2f%%\n", RANGE, (double) (discords_ff_range * 100.0f) / K);
  printf("-----------------------------------------------\n");

  /*--------------------------------------------------------------------------*/
  // RIC imprimo al final del fichero
  int i;
  FILE *fp = fopen(path_result, "a");
  assert(fp);
  //RIC imprimo el tiempo en primer lugar
  fprintf(fp, "#Qsort Time (s): double+float+flexfloat\n");
  fprintf(fp, "%.9f\n", qs_time);
  fprintf(fp, "#Top-100 motifs double (val,i,j)\n");
  for (i = 0; i < K - 1; i++) fprintf(fp, "%f,", sqrt(profile_d[i].val));
  fprintf(fp, "%f\n", sqrt(profile_d[i].val));
  for (i = 0; i < K - 1; i++) fprintf(fp, "%d,", profile_d[i].i);
  fprintf(fp, "%d\n", profile_d[i].i);
  for (i = 0; i < K - 1; i++) fprintf(fp, "%d,", profile_d[i].j);
  fprintf(fp, "%d\n", profile_d[i].j);
  fprintf(fp, "#Top-%d motifs float (val,i,j)\n", K);
  for (i = 0; i < K - 1; i++) fprintf(fp, "%f,", sqrt(profile_f[i].val));
  fprintf(fp, "%f\n", sqrt(profile_f[i].val));
  for (i = 0; i < K - 1; i++) fprintf(fp, "%d,", profile_f[i].i);
  fprintf(fp, "%d\n", profile_f[i].i);
  for (i = 0; i < K - 1; i++) fprintf(fp, "%d,", profile_f[i].j);
  fprintf(fp, "%d\n", profile_f[i].j);
  fprintf(fp, "#Top-100 motifs flexfloat (val,i,j)\n");
  for (i = 0; i < K - 1; i++) fprintf(fp, "%f,", sqrt(ff_get_double(&(profile_ff[i].val))));
  fprintf(fp, "%f\n", sqrt(ff_get_double(&(profile_ff[i].val))));
  for (i = 0; i < K - 1; i++) fprintf(fp, "%d,", profile_ff[i].i);
  fprintf(fp, "%d\n", profile_ff[i].i);
  for (i = 0; i < K - 1; i++) fprintf(fp, "%d,", profile_ff[i].j);
  fprintf(fp, "%d\n", profile_ff[i].j);
  //RIC discords
  fprintf(fp, "#Top-%d discords double (val,i,j)\n", K);
  for (i = ProfileLength - 1; i >= (ProfileLength - K + 1); i--) fprintf(fp, "%f,", sqrt(profile_d[i].val));
  fprintf(fp, "%f\n", sqrt(profile_d[i].val));
  for (i = ProfileLength - 1; i >= (ProfileLength - K + 1); i--) fprintf(fp, "%d,", profile_d[i].i);
  fprintf(fp, "%d\n", profile_d[i].i);
  for (i = ProfileLength - 1; i >= (ProfileLength - K + 1); i--) fprintf(fp, "%d,", profile_d[i].j);
  fprintf(fp, "%d\n", profile_d[i].j);
  fprintf(fp, "#Top-%d discords float (val,i,j)\n", K);
  for (i = ProfileLength - 1; i >= (ProfileLength - K + 1); i--) fprintf(fp, "%f,", sqrt(profile_f[i].val));
  fprintf(fp, "%f\n", sqrt(profile_f[i].val));
  for (i = ProfileLength - 1; i >= (ProfileLength - K + 1); i--) fprintf(fp, "%d,", profile_f[i].i);
  fprintf(fp, "%d\n", profile_f[i].i);
  for (i = ProfileLength - 1; i >= (ProfileLength - K + 1); i--) fprintf(fp, "%d,", profile_f[i].j);
  fprintf(fp, "%d\n", profile_f[i].j);
  fprintf(fp, "#Top-%d discords flexfloat (val,i,j)\n", K);
  for (i = ProfileLength - 1; i >= (ProfileLength - K + 1); i--) fprintf(fp, "%f,", sqrt(ff_get_double(&(profile_ff[i].val))));
  fprintf(fp, "%f\n", sqrt(ff_get_double(&(profile_ff[i].val))));
  for (i = ProfileLength - 1; i >= (ProfileLength - K + 1); i--) fprintf(fp, "%d,", profile_ff[i].i);
  fprintf(fp, "%d\n", profile_ff[i].i);
  for (i = ProfileLength - 1; i >= (ProfileLength - K + 1); i--) fprintf(fp, "%d,", profile_ff[i].j);
  fprintf(fp, "%d\n", profile_ff[i].j);
  fprintf(fp, "#Top-%d accuracy w.r.t double: Motifs float,float\u00B1%d,flexfloat,flexfloat\u00B1%d\n", K, RANGE, RANGE);
  fprintf(fp, "%.2f,%.2f,%.2f,%.2f\n", (double) (motifs_f * 100.0f) / K, (double) (motifs_f_range * 100.0f) / K,
          (double) (motifs_ff * 100.0f) / K, (double) (motifs_ff_range * 100.0f) / K);
  fprintf(fp, "#Top-%d accuracy w.r.t double: Discords float,float\u00B1%d,flexfloat,flexfloat\u00B1%d\n", K, RANGE, RANGE);
  fprintf(fp, "%.2f,%.2f,%.2f,%.2f\n", (double) (discords_f * 100.0f) / K, (double) (discords_f_range * 100.0f) / K,
          (double) (discords_ff * 100.0f) / K, (double) (discords_ff_range * 100.0f) / K);
  fclose(fp);
}

void scrimp_f(float * tSeries, float * AMean, float * ASigma, int ProfileLength,
        int windowSize, int* idx, profItem_f * profile_f, int exclusionZone, int numThreads) {

  /* Private structures initialization -------------------------------- */
  float * profile_tmp = (float *) malloc(sizeof (float)*ProfileLength * numThreads);
  assert(profile_tmp);
  int * profileIndex_tmp = (int *) malloc(sizeof (int)*ProfileLength * numThreads);
  assert(profileIndex_tmp);
  for (int i = 0; i < ProfileLength * numThreads; i++) profile_tmp[i] = INFINITY;
  /* ------------------------------------------------------------------ */

#pragma omp parallel
  {
    float distance, windowSizeDTYPE, lastz;
    int diag, my_offset, i, j, ri;

    windowSizeDTYPE = (float) windowSize;
    my_offset = omp_get_thread_num() * ProfileLength;
    //RIC limpio los flags de cada thread antes de empezar los cálculos
    feclearexcept(FE_ALL_EXCEPT);

#pragma omp for schedule(dynamic)
    for (ri = 0; ri < (ProfileLength - (exclusionZone + 1)); ri++) {
      diag = idx[ri];
      lastz = 0;

      /* Dot product calculation -------------------------- */
      for (j = diag; j < windowSize + diag; j++)
        lastz += tSeries[j] * tSeries[j - diag];

      j = diag;
      i = 0;

      /* Distance calculation ----------------------------- */
      distance = 2 * (windowSizeDTYPE - (lastz - windowSizeDTYPE * AMean[j] * AMean[i]) / (ASigma[j] * ASigma[i]));
      /* -------------------------------------------------- */

      /* Profile update ----------------------------------- */
      if (distance < profile_tmp[my_offset + j]) {
        profile_tmp[my_offset + j] = distance;
        profileIndex_tmp [my_offset + j] = i;
      }
      if (distance < profile_tmp[my_offset + i]) {
        profile_tmp[my_offset + i] = distance;
        profileIndex_tmp [my_offset + i] = j;
      }
      /* -------------------------------------------------- */
      i = 1;
      for (j = diag + 1; j < ProfileLength; j++) {
        /* Dot product update ----------------------- */
        lastz += (tSeries[j + windowSize - 1] * tSeries[i + windowSize - 1]) - (tSeries[j - 1] * tSeries[i - 1]);
        /* ------------------------------------------ */

        /* Distance calculation --------------------- */
        distance = 2 * (windowSizeDTYPE - (lastz - AMean[j] * AMean[i] * windowSizeDTYPE) / (ASigma[j] * ASigma[i]));
        /* ------------------------------------------ */

        /* Profile update --------------------------- */
        if (distance < profile_tmp[my_offset + j]) {
          profile_tmp[my_offset + j] = distance;
          profileIndex_tmp [my_offset + j] = i;
        }

        if (distance < profile_tmp[my_offset + i]) {
          profile_tmp[my_offset + i] = distance;
          profileIndex_tmp[my_offset + i] = j;
        }
        /* ------------------------------------------ */
        i++;
      }

    }
    //RIC el for tiene una barrera implícita a no ser que se ponga nowait
    //#pragma omp barrier

    /* Final profile reduction ---------------------------------- */
    float min_distance;
    int min_index;

#pragma omp for schedule(static)
    for (int colum = 0; colum < ProfileLength; colum++) {
      //RIC para calcular el mínimo, en lugar de poner infinity pongo la primera distancia
      min_distance = INFINITY;

      min_index = 0;

      for (int row = 0; row < numThreads; row++) {
        if (profile_tmp[colum + (row * ProfileLength)] < min_distance) {
          min_distance = profile_tmp[colum + (row * ProfileLength)];
          min_index = profileIndex_tmp[colum + (row * ProfileLength)];
        }
      }
      profile_f[colum].val = min_distance;
      profile_f[colum].i = min_index;
      profile_f[colum].j = colum;
    }
    //RIC el for tiene una barrera implícita a no ser que se ponga nowait
    //#pragma omp barrier
    /* ---------------------------------------------------------- */
    SAVEEXC(omp_get_thread_num());
  }
  free(profile_tmp);
  free(profileIndex_tmp);
}

void scrimp_d(double * tSeries, double * AMean, double * ASigma,
        int ProfileLength, int windowSize, int* idx, profItem_d * profile_d,
        int exclusionZone, int numThreads) {

  /* Private structures initialization -------------------------------- */
  double * profile_tmp = (double *) malloc(sizeof (double) * ProfileLength * numThreads);
  assert(profile_tmp);
  int * profileIndex_tmp = (int *) malloc(sizeof (int) * ProfileLength * numThreads);
  assert(profileIndex_tmp);
  for (int i = 0; i < ProfileLength * numThreads; i++)
    profile_tmp[i] = INFINITY;
  /* ------------------------------------------------------------------ */

#pragma omp parallel
  {
    double distance, windowSizeDTYPE, lastz;
    int diag, my_offset, i, j, ri;

    windowSizeDTYPE = (double) windowSize;
    my_offset = omp_get_thread_num() * ProfileLength;
    //RIC limpio los flags de cada thread antes de empezar los cálculos
    feclearexcept(FE_ALL_EXCEPT);

#pragma omp for schedule(dynamic)
    for (ri = 0; ri < (ProfileLength - (exclusionZone + 1)); ri++) {
      diag = idx[ri];
      lastz = 0;

      /* Dot product calculation -------------------------- */
      for (j = diag; j < windowSize + diag; j++)
        lastz += tSeries[j] * tSeries[j - diag];

      j = diag;
      i = 0;

      /* Distance calculation ----------------------------- */
      distance = 2 * (windowSizeDTYPE - (lastz - windowSizeDTYPE * AMean[j] * AMean[i]) / (ASigma[j] * ASigma[i]));
      /* -------------------------------------------------- */

      /* Profile update ----------------------------------- */
      if (distance < profile_tmp[my_offset + j]) {
        profile_tmp[my_offset + j] = distance;
        profileIndex_tmp [my_offset + j] = i;
      }
      if (distance < profile_tmp[my_offset + i]) {
        profile_tmp[my_offset + i] = distance;
        profileIndex_tmp [my_offset + i] = j;
      }
      /* -------------------------------------------------- */
      i = 1;
      for (j = diag + 1; j < ProfileLength; j++) {
        /* Dot product update ----------------------- */
        lastz += (tSeries[j + windowSize - 1] * tSeries[i + windowSize - 1]) - (tSeries[j - 1] * tSeries[i - 1]);
        /* ------------------------------------------ */

        /* Distance calculation --------------------- */
        distance = 2 * (windowSizeDTYPE - (lastz - AMean[j] * AMean[i] * windowSizeDTYPE) / (ASigma[j] * ASigma[i]));
        /* ------------------------------------------ */

        /* Profile update --------------------------- */
        if (distance < profile_tmp[my_offset + j]) {
          profile_tmp[my_offset + j] = distance;
          profileIndex_tmp [my_offset + j] = i;
        }

        if (distance < profile_tmp[my_offset + i]) {
          profile_tmp[my_offset + i] = distance;
          profileIndex_tmp[my_offset + i] = j;
        }
        /* ------------------------------------------ */
        i++;
      }

    }
    //RIC el for tiene una barrera implícita a no ser que se ponga nowait
    //#pragma omp barrier

    /* Final profile reduction ---------------------------------- */
    double min_distance;
    int min_index;

#pragma omp for schedule(static)
    for (int colum = 0; colum < ProfileLength; colum++) {
      min_distance = INFINITY;

      min_index = 0;

      for (int row = 0; row < numThreads; row++) {
        if (profile_tmp[colum + (row * ProfileLength)] < min_distance) {
          min_distance = profile_tmp[colum + (row * ProfileLength)];
          min_index = profileIndex_tmp[colum + (row * ProfileLength)];
        }
      }
      //RIC meto el sqrt aquí 
      profile_d[colum].val = min_distance;
      profile_d[colum].i = min_index;
      profile_d[colum].j = colum;
    }
    //RIC el for tiene una barrera implícita a no ser que se ponga nowait
    //#pragma omp barrier
    /* ---------------------------------------------------------- */
    SAVEEXC(omp_get_thread_num());
  }
  free(profile_tmp);
  free(profileIndex_tmp);
}

void scrimp_ff(flexfloat_t * tSeries, flexfloat_t* AMean, flexfloat_t* ASigma,
        int ProfileLength, int windowSize, int* idx, profItem_ff * profile_ff,
        int exclusionZone, int numThreads) {

  flexfloat_t windowSize_ff;

  ff_init_int(&windowSize_ff, windowSize, (flexfloat_desc_t) {dist_exp, dist_man});
  /* Private structures initilization --------------------------------- */
  flexfloat_t * profile_priv = (flexfloat_t *) malloc(sizeof (flexfloat_t) * ProfileLength * numThreads);
  assert(profile_priv);
  int * profileIdxs_priv = (int *) malloc(sizeof (int) * ProfileLength * numThreads);
  assert(profileIdxs_priv);

  for (int i = 0; i < ProfileLength * numThreads; i++)
    ff_init_double(&profile_priv[i], INFINITY, (flexfloat_desc_t) {prof_exp, prof_man});
  /* ------------------------------------------------------------------ */

#pragma omp parallel
  {
    flexfloat_t substr, distance, sigma_prods, mean_prods, constant_2;
    flexfloat_t lastz, lastz_cast, dist_cast, mean_cast, sigma_cast;
    uint32_t tid = omp_get_thread_num();
    assert(tid >= 0 && tid < numThreads);
    //printf("mytid: %d\n", tid);
    int myoffset = tid * ProfileLength;

    ff_init_double(&substr, 0.0, (flexfloat_desc_t) {dotp_exp, dotp_man});
    ff_init_double(&distance, 0.0, (flexfloat_desc_t) {dist_exp, dist_man});
    ff_init_double(&mean_prods, 0.0, (flexfloat_desc_t) {stats_exp, stats_man});
    ff_init_double(&sigma_prods, 0.0, (flexfloat_desc_t) {stats_exp, stats_man});
    ff_init_double(&constant_2, 2.0, (flexfloat_desc_t) {dist_exp, dist_man});

    //RIC limpio los flags de cada thread antes de empezar los cálculos
    feclearexcept(FE_ALL_EXCEPT);
#pragma omp for schedule(dynamic)
    for (int ri = 0; ri < (ProfileLength - (exclusionZone + 1)); ri++) {
      int subseq = idx[ri];

      /* Dot product calculation -------------------------- */
      ff_init_double(&lastz, 0.0, (flexfloat_desc_t) {dotp_exp, dotp_man});
      for (int w = 0; w < windowSize; w++)
        ff_fma(&lastz, &tSeries[w + subseq], &tSeries[w], &lastz, tid);

      ff_cast(&lastz_cast, &lastz, (flexfloat_desc_t) {dist_exp, dist_man}, tid);
      /* -------------------------------------------------- */

      /* Distance calculation ----------------------------- */
      ff_mul(&sigma_prods, &ASigma[subseq], &ASigma[0], tid);
      ff_mul(&mean_prods, &AMean[subseq], &AMean[0], tid);

      ff_cast(&mean_cast, &mean_prods, (flexfloat_desc_t) {dist_exp, dist_man}, tid);

      ff_cast(&sigma_cast, &sigma_prods, (flexfloat_desc_t) {dist_exp, dist_man}, tid);
      ff_mul(&distance, &mean_cast, &windowSize_ff, tid);
      ff_sub(&distance, &lastz_cast, &distance, tid);
      ff_div(&distance, &distance, &sigma_cast, tid);
      ff_sub(&distance, &windowSize_ff, &distance, tid);
      ff_mul(&distance, &distance, &constant_2, tid);
      //RIC meto el ABS
      ABS_FF(&distance, ((flexfloat_desc_t) {dist_exp, dist_man}));
      //ZERO_FF(&distance, ((flexfloat_desc_t) {dist_exp, dist_man}));
      /* -------------------------------------------------- */

      /* Profile update ----------------------------------- */
      ff_cast(&dist_cast, &distance, (flexfloat_desc_t) {prof_exp, prof_man}, tid);
      if (ff_lt(&dist_cast, &profile_priv[subseq + myoffset], tid)) {
        profile_priv[subseq + myoffset] = dist_cast;
        profileIdxs_priv[subseq + myoffset] = 0;
      }
      if (ff_lt(&dist_cast, &profile_priv[myoffset], tid)) {
        profile_priv[myoffset] = dist_cast;
        profileIdxs_priv[myoffset] = subseq;
      }
      /* -------------------------------------------------- */
      int i = 1;

      for (int j = subseq + 1; j < ProfileLength; j++) {
        /* Dot product update ----------------------- */
        ff_fma(&lastz, &tSeries[j + windowSize - 1], &tSeries[i + windowSize - 1], &lastz, tid);
        ff_mul(&substr, &tSeries[j - 1], &tSeries[i - 1], tid);
        ff_sub(&lastz, &lastz, &substr, tid);
        ff_cast(&lastz_cast, &lastz, (flexfloat_desc_t) {dist_exp, dist_man}, tid);
        /* ------------------------------------------ */

        /* Distance calculation --------------------- */
        ff_mul(&sigma_prods, &ASigma[j], &ASigma[i], tid);
        ff_mul(&mean_prods, &AMean[j], &AMean [i], tid);
        ff_cast(&mean_cast, &mean_prods, (flexfloat_desc_t) {dist_exp, dist_man}, tid);
        ff_cast(&sigma_cast, &sigma_prods, (flexfloat_desc_t) {dist_exp, dist_man}, tid);
        ff_mul(&distance, &mean_cast, &windowSize_ff, tid);
        ff_sub(&distance, &lastz_cast, &distance, tid);
        ff_div(&distance, &distance, &sigma_cast, tid);
        ff_sub(&distance, &windowSize_ff, &distance, tid);
        ff_mul(&distance, &distance, &constant_2, tid);
        //RIC meto el ABS
        ABS_FF(&distance, ((flexfloat_desc_t) {dist_exp, dist_man}));
        //ZERO_FF(&distance,((flexfloat_desc_t) {dist_exp, dist_man}));
        /* ------------------------------------------ */

        /* Profile update --------------------------- */
        ff_cast(&dist_cast, &distance, (flexfloat_desc_t) {prof_exp, prof_man}, tid);
        if (ff_lt(&dist_cast, &profile_priv[j + myoffset], tid)) {
          profile_priv[j + myoffset] = dist_cast;
          profileIdxs_priv[j + myoffset] = i;
        }
        if (ff_lt(&dist_cast, &profile_priv[i + myoffset], tid)) {
          profile_priv[i + myoffset] = dist_cast;
          profileIdxs_priv[i + myoffset] = j;
        }
        /* ------------------------------------------ */
        i++;
      }

    }
    //RIC el for tiene una barrera implícita a no ser que se ponga nowait
    //#pragma omp barrier

    /* Final profile reduction ---------------------------------- */
    flexfloat_t min_distance;
    int min_index;
    fexcept_t flagsPriv;
#pragma omp for schedule(static)
    for (int colum = 0; colum < ProfileLength; colum++) {
      //RIC parece que la asignación a INFINITY levanta el flag de overflow en flexfloat
      fegetexceptflag(&flagsPriv, FE_ALL_EXCEPT);

      ff_init_double(&min_distance, INFINITY, (flexfloat_desc_t) {prof_exp, prof_man});
      fesetexceptflag(&flagsPriv, FE_ALL_EXCEPT);

      min_index = 0;

      for (int row = 0; row < numThreads; row++) {
        if (ff_lt(&profile_priv[colum + (row * ProfileLength)], &min_distance, tid)) {
          min_distance = profile_priv[colum + (row * ProfileLength)];
          min_index = profileIdxs_priv[colum + (row * ProfileLength)];
        }
      }
      profile_ff[colum].val = min_distance;
      profile_ff[colum].i = min_index;
      profile_ff[colum].j = colum;
    }
    //RIC el for tiene una barrera implícita a no ser que se ponga nowait
    //#pragma omp barrier
    /* ---------------------------------------------------------- */
    SAVEEXC(omp_get_thread_num());
  }
  free(profile_priv);
  free(profileIdxs_priv);
}

void read_config(char * file_name) {

  char * path_config = malloc(1000 * sizeof (char));

  strcpy(path_config, PATH_CFG);
  strcat(path_config, file_name);
  path_config[strlen(path_config) - 4] = 0;
  strcat(path_config, ".cfg");

  FILE * file = fopen(path_config, "r");
  int err1, err2, err3, err4, err5, err6, err7, err8;

  if (file == NULL) printf("CFG FILE ERRROR\n");
  err1 = fscanf(file, "%d", &dist_exp);
  err2 = fscanf(file, "%d", &dist_man);
  err3 = fscanf(file, "%d", &dotp_exp);
  err4 = fscanf(file, "%d", &dotp_man);
  err5 = fscanf(file, "%d", &stats_exp);
  err6 = fscanf(file, "%d", &stats_man);
  err7 = fscanf(file, "%d", &prof_exp);
  err8 = fscanf(file, "%d", &prof_man);

  if (!err1 || !err2 || !err3 || !err4 || !err5
          || !err6 || !err7 || !err8)
    printf("CFG FILE ERRROR\n");

  free(path_config);
  fclose(file);
}

void print_header() {
  /* Printing program header ------------------------------------------ */
  printf("\n");
  printf( "███████╗ ██████╗██████╗ ██╗███╗   ███╗██████╗ \n"
          "██╔════╝██╔════╝██╔══██╗██║████╗ ████║██╔══██╗\n"
          "███████╗██║     ██████╔╝██║██╔████╔██║██████╔╝\n"
          "╚════██║██║     ██╔══██╗██║██║╚██╔╝██║██╔═══╝ \n"
          "███████║╚██████╗██║  ██║██║██║ ╚═╝ ██║██║     \n"
          "╚══════╝ ╚═════╝╚═╝  ╚═╝╚═╝╚═╝     ╚═╝╚═╝     \n"
          "+-+-+-+-+-+-+-+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-+\n"
          "| F l e x F l o a t        B e n c h m a r k |\n"
          "+-+-+-+-+-+-+-+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-+\n");
  printf("\n");
  /* ------------------------------------------------------------------ */
}

