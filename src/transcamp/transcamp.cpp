#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <vector>
#include <string>
#include <sstream>
#include <chrono>
#include <algorithm> //RIC para sort()
#include <assert.h>
#include <omp.h>
#include <unistd.h> //RIC para getpid(), usado para obtener un nombre de fichero único
#include "flexfloat.hpp" //RIC para la librería flexfloat
extern "C"
{
#include "rapl-read.h"
}


//RIC para no tener que poner std:: todo el rato
using namespace std;

#define PATH_CFG "./configs/"
#define PATH_TIME_SERIES "../../tseries/"
#define PATH_RESULTS "./"


int numThreads, exclusionZone, windowSize, timeSeriesLength, profileLength;
uint8_t covExp, covMan, corrExp, corrMan, statsExp, statsMan, profExp, profMan;
double scaleFactor;

// Computes all required statistics for SCAMP, populating info with these values

void preprocess(vector<double> &tSeries_d, vector<flexfloat_t> &means_ff, vector<flexfloat_t> &norms_ff,
        vector<flexfloat_t> &df_ff, vector<flexfloat_t> &dg_ff) {
  vector<double> prefix_sum(tSeries_d.size());
  vector<double> prefix_sum_sq(tSeries_d.size());
  vector<double> means_d(tSeries_d.size());
  vector<double> norms_d(tSeries_d.size());
  vector<double> df_d(tSeries_d.size());
  vector<double> dg_d(tSeries_d.size());  

//RIC esto es el profileLength que está definido como variable global
  //int n = tSeries_d.size() - windowSize + 1;

  //RIC se calculan los vectores de suma prefija y cuadrado de la suma
  prefix_sum [0] = tSeries_d[0];
  prefix_sum_sq [0] = tSeries_d[0] * tSeries_d[0];
  for (int i = 1; i < tSeries_d.size(); ++i) {
    prefix_sum[i] = tSeries_d[i] + prefix_sum[i - 1];
    prefix_sum_sq[i] = tSeries_d[i] * tSeries_d[i] + prefix_sum_sq[i - 1];
  }

  //RIC se utilizan las sumas para calcular la media de los valores de cada ventana como
  // el último valor de la ventana menos el primero partido por el tamaño de la ventana
  //CHECK_STATS(means_d[0]);
  //RIC instancio los vectores float y flexfloat

  ff_init_double(&means_ff[0], means_d[0], (flexfloat_desc_t) {statsExp, statsMan});
  for (int i = 1; i < profileLength; ++i) {
    //CHECK_STATS(means_d[i]);
    //RIC instancio los vectores float y flexfloat
    ff_init_double(&means_ff[i], means_d[i], (flexfloat_desc_t) {statsExp, statsMan});
  }

  double sum = 0;
  for (int i = 0; i < windowSize; ++i) {
    double val = tSeries_d[i] - means_d[0];
    sum += val * val;
  }
  norms_d[0] = sum;
  
  //RIC se calculan las L2-norms (norma euclídea, distancia euclídea)
  for (int i = 1; i < profileLength; ++i)
    norms_d[i] = norms_d[i - 1] + ((tSeries_d[i - 1] - means_d[i - 1]) + (tSeries_d[i + windowSize - 1] - means_d[i])) * 
            (tSeries_d[i + windowSize - 1] - tSeries_d[i - 1]);
  //RIC The default type for floating-point literals is double (http://www.cplusplus.com/doc/tutorial/constants/)
  // F para floats
  for (int i = 0; i < profileLength; ++i) {
    norms_d[i] = 1.0 / sqrt(norms_d[i]);
    ff_init_double(&norms_ff[i], norms_d[i], (flexfloat_desc_t) {statsExp, statsMan});
  }

  //RIC Comprobar esto!! Por qué profileLength-1 si luego se usa hasta profileLength en scamp
  for (int i = 0; i < profileLength - 1; ++i) {
    df_d[i] = (tSeries_d[i + windowSize] - tSeries_d[i]) / 2.0;
    dg_d[i] = (tSeries_d[i + windowSize] - means_d[i + 1]) + (tSeries_d[i] - means_d[i]);

    ff_init_double(&df_ff[i], df_d[i], (flexfloat_desc_t) {statsExp, statsMan});
    ff_init_double(&dg_ff[i], dg_d[i], (flexfloat_desc_t) {statsExp, statsMan});
  }
}


void scamp_ff(vector<flexfloat_t> &tSeries, vector<flexfloat_t> &means, vector<flexfloat_t> &norms,
        vector<flexfloat_t> &df, vector<flexfloat_t> &dg, vector<flexfloat_t> &profile, vector<unsigned> profileIndex){
  /* Private structures initialization -------------------------------- */
  vector<flexfloat_t> profile_tmp(profileLength * numThreads);
  vector<int> profileIndex_tmp(profileLength * numThreads);

  for (int i = 0; i < profileLength * numThreads; i++)
    ff_init_double(&profile_tmp[i], -numeric_limits<double>::infinity(), (flexfloat_desc_t) {profExp, profMan});
  /* ------------------------------------------------------------------ */

#pragma omp parallel
  {
    uint32_t tid = omp_get_thread_num();
    assert(tid >= 0 && tid < numThreads);
    int my_offset = tid * profileLength;
    flexfloat_t covCov, tmpCov, tmp1Cov, tmp2Cov, covCorr, corrCorr, tmp1Corr, tmp2Corr, corrProf;

    ff_init_double(&tmpCov, 0.0, (flexfloat_desc_t) {covExp, covMan});
    ff_init_double(&tmp1Cov, 0.0, (flexfloat_desc_t) {covExp, covMan});
    ff_init_double(&tmp2Cov, 0.0, (flexfloat_desc_t) {covExp, covMan});
    ff_init_double(&corrCorr, 0.0, (flexfloat_desc_t) {corrExp, corrMan});
    ff_init_double(&covCorr, 0.0, (flexfloat_desc_t) {corrExp, corrMan});
    ff_init_double(&tmp1Corr, 0.0, (flexfloat_desc_t) {corrExp, corrMan});
    ff_init_double(&tmp2Corr, 0.0, (flexfloat_desc_t) {corrExp, corrMan});
    
#pragma omp for schedule(dynamic)
    // Go through diagonals
    for (int diag = exclusionZone + 1; diag < profileLength; diag++) {

      ff_init_double(&covCov, 0.0, (flexfloat_desc_t) {covExp, covMan});
      for (int i = 0; i < windowSize; i++) {
        //covariance += ((tSeries[diag + i] - means[diag]) * (tSeries[i] - means[0]));
        ff_cast(&tmp1Cov, &means[diag], (flexfloat_desc_t) {covExp, covMan}, tid);
        ff_cast(&tmp2Cov, &means[0], (flexfloat_desc_t) {covExp, covMan}, tid);
        ff_sub(&tmp1Cov, &tSeries[diag + i], &tmp1Cov, tid);
        ff_sub(&tmp2Cov, &tSeries[i], &tmp2Cov, tid);
        //Cambio la multiplicación y la suma por un fma, operación que existe en las FPU (Floating Point Unit)
        //ff_mul(&tmp, &tmp1, &tmp2, tid);
        //ff_acc(&covariance, &tmp, tid);
        //RIC ff_fma(flexfloat_t *dest, const flexfloat_t *a, const flexfloat_t *b, const flexfloat_t *c, const uint32_t tid)
        // multiplica a*b le suma c y lo deja en dest
        ff_fma(&covCov, &tmp1Cov, &tmp2Cov, &covCov, tid);
      }

      int i = 0;
      int j = diag;
      
      //RIC pruebo poniendo la parte del cálculo gordo del principio a high y los updates a low
      ff_cast(&covCorr, &covCov, (flexfloat_desc_t) {corrExp, corrMan}, tid);
      //correlation = covariance * norms[offset] * norms[j];
      ff_cast(&tmp1Corr, &norms[i], (flexfloat_desc_t) {corrExp, corrMan}, tid);
      ff_cast(&tmp2Corr, &norms[j], (flexfloat_desc_t) {corrExp, corrMan}, tid);
      ff_mul(&corrCorr, &covCorr, &tmp1Corr, tid);
      ff_mul(&corrCorr, &corrCorr, &tmp2Corr, tid);
      //la asociatividad puede que no de el mismo resultado (mejor multiplicar primero por una norma y luego por la otra)
      //ff_mul(&correLo, &norms[i], &norms[j], tid);
      //ff_mul(&correLo, &correLo, &covLo, tid);
      
      ff_cast(&corrProf, &corrCorr, (flexfloat_desc_t) {profExp, profMan}, tid);
      //Profile update
      
      //if (correlation > profile_tmp[offset + my_offset]) {
      if (ff_gt(&corrProf, &profile_tmp[i + my_offset], tid)) {
        profile_tmp[i + my_offset] = corrProf;
        profileIndex_tmp[i + my_offset] = j;
      }
      //if (correlation > profile_tmp[j + my_offset]) {
      if (ff_gt(&corrProf, &profile_tmp[j + my_offset], tid)) {
        profile_tmp[j + my_offset] = corrProf;
        profileIndex_tmp[j + my_offset] = i;
      }
      i = 1;
      for (int j = diag + 1; j < profileLength; j++) {
        //Distancia?
        //covariance += (df[offset-1] * dg[j-1] + df[j-1] * dg[offset-1]);
        ff_cast(&tmp1Cov, &df[i-1], (flexfloat_desc_t) {covExp, covMan}, tid);
        ff_cast(&tmp2Cov, &dg[j-1], (flexfloat_desc_t) {covExp, covMan}, tid);
        ff_mul(&tmpCov, &tmp1Cov, &tmp2Cov, tid);
        ff_cast(&tmp1Cov, &df[j-1], (flexfloat_desc_t) {covExp, covMan}, tid);
        ff_cast(&tmp2Cov, &dg[i-1], (flexfloat_desc_t) {covExp, covMan}, tid);
        ff_mul(&tmp1Cov, &tmp1Cov, &tmp2Cov, tid);
        ff_add(&tmpCov, &tmp1Cov, &tmpCov, tid);
        ff_acc(&covCov, &tmpCov, tid);
        
        ff_cast(&covCorr, &covCov, (flexfloat_desc_t) {corrExp, corrMan}, tid);
        //correlation = covariance * norms[offset] * norms[j];
        ff_cast(&tmp1Corr, &norms[i], (flexfloat_desc_t) {corrExp, corrMan}, tid);
        ff_cast(&tmp2Corr, &norms[j], (flexfloat_desc_t) {corrExp, corrMan}, tid);
        ff_mul(&corrCorr, &covCorr, &tmp1Corr, tid);
        ff_mul(&corrCorr, &corrCorr, &tmp2Corr, tid);
        
        ff_cast(&corrProf, &corrCorr, (flexfloat_desc_t) {profExp, profMan}, tid);
        //Profile update
        //if (correlation > profile_tmp[offset + my_offset]) {
        if (ff_gt(&corrProf, &profile_tmp[i + my_offset], tid)) {
          profile_tmp[i + my_offset] = corrProf;
          profileIndex_tmp[i + my_offset] = j;
        }
        //if (correlation > profile_tmp[j + my_offset]) {
        if (ff_gt(&corrProf, &profile_tmp[j + my_offset], tid)) {
          profile_tmp[j + my_offset] = corrProf;
          profileIndex_tmp[j + my_offset] = i;
        }
        i++;
      }
    }

    flexfloat_t max_corr;
    int max_index;

#pragma omp for schedule(static)
    for (int colum = 0; colum < profileLength; colum++) {
      //max_corr = -numeric_limits<flexfloat_t>::infinity();

      ff_init_double(&max_corr, -numeric_limits<double>::infinity(), (flexfloat_desc_t) {profExp, profMan});
      for (int row = 0; row < numThreads; row++) {
        //if (profile_tmp[colum + (row * profileLength)] > max_corr) {
        if (ff_gt(&profile_tmp[colum + (row * profileLength)], &max_corr, tid)) {
          max_corr = profile_tmp[colum + (row * profileLength)];
          max_index = profileIndex_tmp[colum + (row * profileLength)];
        }
      }
      profile[colum] = max_corr;
      profileIndex[colum] = max_index;
    }
    //RIC el for tiene una barrera implícita a no ser que se ponga nowait
    //#pragma omp barrier
  }
}

  /*string pathcfg(PATH_CFG);
  pathcfg += filename.substr(0, filename.size() - 4) + ".cfg";

  fstream cfgFile(pathcfg, ios_base::in);
  if (!cfgFile.is_open()) {
    cout << "No se pudo abrir el archivo de configuración flexfloat: " << pathcfg << endl;
    exit(EXIT_FAILURE);
  }
  int tmp;
  cfgFile >> tmp;
  covExp = tmp;
  cfgFile >> tmp;
  covMan = tmp;
  cfgFile >> tmp;
  corrExp = tmp;
  cfgFile >> tmp;
  corrMan = tmp;
  cfgFile >> tmp;
  statsExp = tmp;
  cfgFile >> tmp;
  statsMan = tmp;
  cfgFile >> tmp;
  profExp = tmp;
  cfgFile >> tmp;
  profMan = tmp;

  if (!covExp || !covMan || !profExp || !profMan || !corrExp || !corrMan || !statsExp || !statsMan) {
    cout << "Error en el archivo de configuración. Algún valor es 0." << endl;
    exit(EXIT_FAILURE);
  }
  cfgFile.close();
}*/

int main(int argc, char* argv[]) {
  try {
    // Creation of time meassure structures
    chrono::high_resolution_clock::time_point tstart, tend;
    chrono::duration<double> time_elapsed, mp_time_f, mp_time_d, mp_time_ff;

    if (argc != 7) {
      cout << "[ERROR] usage: ./scamp tseries.txt win_size num_threads scale_factor exponent mantissa" << endl;
      //cout << "[ERROR] usage2: ./scamp tseries.txt win_size num_threads scale_factor error_diffusion FF_exponent FF_mantissa" << endl;
      return 0;
    }

    windowSize = atoi(argv[2]);
    numThreads = atoi(argv[3]);
    scaleFactor = atof(argv[4]);
    //if (argc == 6) {
    //  readConfig(argv[1]);
    //} else {
      covExp = atoi(argv[5]);
      covMan = atoi(argv[6]);
      corrExp = covExp;
      corrMan = covMan;
      statsExp = covExp;
      statsMan = covMan;
      profExp = covExp;
      profMan = covMan;
    //}
    omp_set_num_threads(numThreads);

    // Set the exclusion zone
    exclusionZone = (int) (windowSize * 0.25);

    //RIC creo los vectores que almacenarán la serie con sus distintos tipos
    //    una vez sé el exponente y la mantisa para flexfloat
    vector<float> tSeries_f;
    vector<double> tSeries_d;
    vector<flexfloat_t> tSeries_ff;

    ff_start_stats(numThreads);

    string inputfilename = argv[1];
    //RIC meto en un stringstream el nombre del fichero de salida con todos los parámetros
    stringstream tmp;
    tmp << PATH_RESULTS << inputfilename.substr(0, inputfilename.size() - 4) << "_w" <<
            windowSize << "_s" << setprecision(2) << scaleFactor <<  "_t" << numThreads <<
            "_cov" << (int) covExp << "-" << (int) covMan << "_corr" << (int) corrExp << "-" << (int) corrMan << 
            "_stat" << (int) statsExp << "-" << (int) statsMan << "_prof" << (int) profExp << "-" << (int) profMan <<
            "_" << getpid() << ".csv";
    string outfilename = tmp.str();

    rapl_perf_init();

    // Display info through console
    cout << endl;
    cout << "############################################################" << endl;
    cout << "/////////////////////// TranSCAMP //////////////////////////" << endl;
    cout << "############################################################" << endl;
    cout << endl;
    cout << "[>>] Reading File: " << inputfilename << "..." << endl;

    /* ------------------------------------------------------------------ */
    /* Read time series file */
    tstart = chrono::high_resolution_clock::now();

    fstream timeSeriesFile(string(PATH_TIME_SERIES) + inputfilename, ios_base::in);
    double tempval, lastQuantError = 0.0, tSeriesMin = numeric_limits<double>::infinity(), tSeriesMax = -numeric_limits<double>::infinity();
    timeSeriesLength = 0;
    while (timeSeriesFile >> tempval) {
      tempval *= scaleFactor;
      tSeries_d.push_back(tempval);
      flexfloat_t tmp;
      ff_init_double(&tmp, tempval, (flexfloat_desc_t) {covExp, covMan});
      tSeries_ff.push_back(tmp);
      
      timeSeriesLength++;
    }
    timeSeriesFile.close();

    tend = chrono::high_resolution_clock::now();
    time_elapsed = tend - tstart;
    cout << "[OK] Read File Time: " << setprecision(numeric_limits<double>::digits10 + 2) << time_elapsed.count() << " seconds." << endl;

    // Set Matrix Profile Length
    profileLength = timeSeriesLength - windowSize + 1;

    //RIC Una vez leída la serie ya tengo todos los parámetros necesarios para crear los vectores auxiliares
    vector<flexfloat_t> norms_ff(timeSeriesLength), means_ff(timeSeriesLength), df_ff(timeSeriesLength), dg_ff(timeSeriesLength);
    vector<flexfloat_t> profile_ff(profileLength);
    vector<unsigned> profileIndex(profileLength);

    // Display info through console
    cout << endl;
    cout << "------------------------------------------------------------" << endl;
    cout << "************************** INFO ****************************" << endl;
    cout << endl;
    cout << " Time series length:  " << timeSeriesLength << endl;
    cout << " Window size:         " << windowSize << endl;
    cout << " Scale factor:        " << setprecision(4) << scaleFactor << endl;
    cout << " Number of threads:   " << numThreads << endl;
    cout << " Exclusion zone:      " << exclusionZone << endl;
    cout << " Profile length:      " << profileLength << endl;
    cout << " FF cov  prec exp, man: " << (int) covExp << ", " << (int) covMan << endl;
    cout << " FF corr prec exp, man: " << (int) corrExp << ", " << (int) corrMan << endl;
    cout << " FF stat prec exp, man: " << (int) statsExp << ", " << (int) statsMan << endl;
    cout << " FF prof prec exp, man: " << (int) profExp << ", " << (int) profMan << endl;
    cout << "------------------------------------------------------------" << endl;
    cout << endl;

    // Preprocess, statistics, get the mean and standard deviation of every subsequence in the time series
    cout << "[>>] Preprocessing..." << endl;

    tstart = chrono::high_resolution_clock::now();
    preprocess(tSeries_d, means_ff, norms_ff, df_ff, dg_ff);
    tend = chrono::high_resolution_clock::now();
    time_elapsed = tend - tstart;
    cout << "[OK] Preprocess Time:         " << setprecision(numeric_limits<double>::digits10 + 2) <<
            time_elapsed.count() << " seconds." << endl;

    /***********************************************/
    /******************** SCAMP flexfloat ********************/
    cout << "[>>] Performing SCAMP flexfloat..." << endl;
    tstart = chrono::high_resolution_clock::now();
    rapl_perf_start(0);
    scamp_ff(tSeries_ff, means_ff, norms_ff, df_ff, dg_ff, profile_ff, profileIndex);
    rapl_perf_end();
    tend = chrono::high_resolution_clock::now();
    mp_time_ff = tend - tstart;
    cout << "[OK] SCAMP flexfloat Time:              " << setprecision(numeric_limits<float>::digits10 + 2) <<
            mp_time_ff.count() << " seconds." << endl;
    /***********************************************/

    //RIC imprimo el tiempo en primer lugar
    cout << "[>>] Saving Stats: " << outfilename << " ..." << endl;
    fstream statsFile(outfilename, ios_base::out);

    statsFile << "#Time (s): double,float,flexfloat" << endl;
    statsFile << setprecision(9) << mp_time_d.count() << ", " << mp_time_f.count() << ", " << mp_time_ff.count() << endl;
    statsFile << "#Profile Length" << endl;
    statsFile << profileLength << endl;
    statsFile << "#i,tseries,doubleMP,j,floatMP,j,errorf(%),flexfloatMP,j,errorff(%)" << endl;
    for (int i = 0; i < profileLength; i++) {
      statsFile << sqrt(2 * windowSize * (1 - ff_get_double(&(profile_ff[i])))) << " " << profileIndex[i] << endl;
    }
    statsFile.close();
    cout << "------------------------------------------------------------" << endl;
    cout << "********************* FlexFloat stats: *********************" << endl;
    ff_print_stats();
    //ffStatsToFile(outfilename, numThreads);
    cout << "------------------------------------------------------------" << endl;
    /* ------------------------------------------------------------------------ */

    cout << endl;

    rapl_perf_print();

  } catch (exception &e) {
    cout << "Something went wrong: " << e.what() << endl;
  }
}
