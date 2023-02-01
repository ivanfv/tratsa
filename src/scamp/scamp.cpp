#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <vector>
#include <string>
#include <sstream>
#include <chrono>
#include <algorithm>
#include <assert.h>
#include <omp.h>
#include <unistd.h>
#include <typeinfo>

#define PATH_TIME_SERIES "./"
#define PATH_RESULTS "./results/MP_"

#define DTYPE double	/* DATA TYPE */
#define ITYPE unsigned	/* INDEX TYPE */

using namespace std;

int numThreads, exclusionZone, windowSize, timeSeriesLength, profileLength;

// Computes all required statistics for SCAMP, populating info with these values
void preprocess(vector<double> &tSeries, vector<double> &means, vector<double> &norms,
        vector<DTYPE> &df, vector<DTYPE> &dg) {

  vector<double> prefix_sum(tSeries.size());
  vector<double> prefix_sum_sq(tSeries.size());

  // Calculates prefix sum and square sum vectors 
  prefix_sum [0] = tSeries[0];
  prefix_sum_sq [0] = tSeries[0] * tSeries[0];
  for (ITYPE i = 1; i < tSeries.size(); ++i) 
  {
    prefix_sum[i]    = tSeries[i] + prefix_sum[i - 1];
    prefix_sum_sq[i] = tSeries[i] * tSeries[i] + prefix_sum_sq[i - 1];
  }

  // Prefix sum value is used to calculate mean value of a given window, taking last value
  // of the window minus the first one and dividing by window size	
  means[0] = prefix_sum[windowSize - 1] / static_cast<DTYPE> (windowSize);
  for (ITYPE i = 1; i < profileLength; ++i) 
  {
    means[i] = (prefix_sum[i + windowSize - 1] - prefix_sum[i - 1]) / static_cast<double> (windowSize);
  }

  double sum = 0;
  for (ITYPE i = 0; i < windowSize; ++i) 
  {
    double val = tSeries[i] - means[0];
    sum += val * val;
  }
  norms[0] = sum;

  // Calculates L2-norms (euclidean norm, euclidean distance)
  for (ITYPE i = 1; i < profileLength; ++i)
  {
    norms[i] = norms[i - 1] + ((tSeries[i - 1] - means[i - 1]) + (tSeries[i + windowSize - 1] - means[i])) * 
            (tSeries[i + windowSize - 1] - tSeries[i - 1]);
  }
  for (ITYPE i = 0; i < profileLength; ++i) 
  {
    norms[i] = (1.0 / sqrt(norms[i]));
  }

  // Calculates df and dg vectors
  for (ITYPE i = 0; i < profileLength - 1; ++i) {
    df[i] = (DTYPE)(tSeries[i + windowSize] - tSeries[i]) / 2.0;
    dg[i] = (DTYPE)(tSeries[i + windowSize] - means[i + 1]) + (tSeries[i] - means[i]);
  }
}


void scamp(vector<double> &tSeries, vector<double> &means, vector<double> &norms,
        vector<DTYPE> &df, vector<DTYPE> &dg, vector<DTYPE> &profile, vector<ITYPE> &profileIndex) {
  /* Private structures initialization -------------------------------- */
  vector<DTYPE> profile_tmp(profileLength * numThreads);
  vector<int> profileIndex_tmp(profileLength * numThreads);
  for (int i = 0; i < profileLength * numThreads; i++) profile_tmp[i] = -numeric_limits<DTYPE>::infinity();
  /* ------------------------------------------------------------------ */

  #pragma omp parallel
  {
    int my_offset = omp_get_thread_num() * profileLength;
    DTYPE covariance, correlation;
    
    #pragma omp for schedule(dynamic)
    // Go through diagonals
    for (int diag = exclusionZone + 1; diag < profileLength; diag++) 
    {
      covariance = 0;
      for (int i = 0; i < windowSize; i++) 
      {
        covariance += (((DTYPE)tSeries[diag + i] - (DTYPE)means[diag]) * ((DTYPE)tSeries[i] - (DTYPE)means[0]));
      }

      int i = 0;
      int j = diag;

      correlation = covariance * (DTYPE)norms[i] * (DTYPE)norms[j];

      if (correlation > profile_tmp[i + my_offset]) 
      {
        profile_tmp[i + my_offset] = correlation;
        profileIndex_tmp[i + my_offset] = j;
      }
      if (correlation > profile_tmp[j + my_offset]) 
      {
        profile_tmp[j + my_offset] = correlation;
        profileIndex_tmp[j + my_offset] = i;
      }

      i = 1;

      for (int j = diag + 1; j < profileLength; j++) 
      {
		if(i % 65536 == 0)
		{
			covariance = 0;
			for (int w = 0; w < windowSize; w++) 
			{
				covariance += (((DTYPE)tSeries[j + w] - (DTYPE)means[j]) * ((DTYPE)tSeries[w + i] - (DTYPE)means[i]));
			}			
		}
		else
		{
			covariance += (df[i - 1] * dg[j - 1] + df[j - 1] * dg[i - 1]);
		}
		

		
        correlation = covariance * (DTYPE)norms[i] * (DTYPE)norms[j];


        if (correlation > profile_tmp[i + my_offset]) 
        {
          profile_tmp[i + my_offset] = correlation;
          profileIndex_tmp[i + my_offset] = j;
        }

        if (correlation > profile_tmp[j + my_offset]) 
        {
          profile_tmp[j + my_offset] = correlation;
          profileIndex_tmp[j + my_offset] = i;
        }
        i++;
      }
    }
    
    

    DTYPE max_corr;
    ITYPE max_index;

    #pragma omp for schedule(static)
    for (int colum = 0; colum < profileLength; colum++) {
      max_corr = -numeric_limits<DTYPE>::infinity();
      for (int row = 0; row < numThreads; row++) {
        if (profile_tmp[colum + (row * profileLength)] > max_corr) {
          max_corr = profile_tmp[colum + (row * profileLength)];
          max_index = profileIndex_tmp[colum + (row * profileLength)];
        }
      }
      profile[colum] = max_corr;
      profileIndex[colum] = max_index;
    }
  }
}



int main(int argc, char* argv[]) {
  try {
    // Creation of time meassure structures
    chrono::high_resolution_clock::time_point tstart, tend;
    chrono::duration<double> time_elapsed, mp_time_f, mp_time_d, mp_time_ff;

    if (argc != 4) {
      cout << "[ERROR] usage: ./scamp input_file win_size num_threads" << endl;
      return 0;
    }

    windowSize = atoi(argv[2]);
    numThreads = atoi(argv[3]);

    omp_set_num_threads(numThreads);

    // Set the exclusion zone to 0.25
    exclusionZone = (ITYPE) (windowSize * 0.25);

    vector<double> tSeries;

    string inputfilename = argv[1];

    stringstream tmp;
    tmp <<  PATH_RESULTS << inputfilename.substr(0, inputfilename.size() - 4) << "_w" <<
            windowSize << "_t" << numThreads << "_pid" << getpid() << ".txt";
    string outfilename = tmp.str();

    // Display info through console
    cout << endl;
    cout << "############################################################" << endl;
    cout << "///////////////////////// SCAMP ////////////////////////////" << endl;
    cout << "############################################################" << endl;
    cout << endl;
    cout << "[>>] Reading File: " << inputfilename << "..." << endl;

    /* ------------------------------------------------------------------ */
    /* Read time series file */
    tstart = chrono::high_resolution_clock::now();

    fstream timeSeriesFile(PATH_TIME_SERIES + inputfilename, ios_base::in);

    double tempval, tSeriesMin = numeric_limits<DTYPE>::infinity(), tSeriesMax = -numeric_limits<DTYPE>::infinity();
  
    timeSeriesLength = 0;
    while (timeSeriesFile >> tempval) {
      tSeries.push_back(tempval);

      if (tempval < tSeriesMin) tSeriesMin = tempval;
      if (tempval > tSeriesMax) tSeriesMax = tempval;
      timeSeriesLength++;
    }
    timeSeriesFile.close();

//    for(int i = 0; i < timeSeriesLength; i++) tSeries[i] = (tSeries[i] - tSeriesMin) / (tSeriesMax - tSeriesMin);

    tend = chrono::high_resolution_clock::now();
    time_elapsed = tend - tstart;
    cout << "[OK] Read File Time: " << setprecision(numeric_limits<DTYPE>::digits10 + 2) << time_elapsed.count() << " seconds." << endl;

    // Set Matrix Profile Length
    profileLength = timeSeriesLength - windowSize + 1;

    // Auxiliary vectors
    vector<DTYPE>  df(timeSeriesLength), dg(timeSeriesLength);
    vector<double> norms(timeSeriesLength), means(timeSeriesLength);
    vector<DTYPE> profile(profileLength);
    vector<ITYPE> profileIndex(profileLength);


    // Display info through console
    cout << endl;
    cout << "------------------------------------------------------------" << endl;
    cout << "************************** INFO ****************************" << endl;
    cout << endl;
    cout << " Data type:           " << typeid(df).name() << endl;
    cout << " Time series length:  " << timeSeriesLength << endl;
    cout << " Window size:         " << windowSize << endl;
    cout << " Time series min:     " << tSeriesMin << endl;
    cout << " Time series max:     " << tSeriesMax << endl;
    cout << " Number of threads:   " << numThreads << endl;
    cout << " Exclusion zone:      " << exclusionZone << endl;
    cout << " Profile length:      " << timeSeriesLength << endl;
    cout << "------------------------------------------------------------" << endl;
    cout << endl;

    /***************** Preprocess ******************/
    cout << "[>>] Preprocessing..." << endl;
    tstart = chrono::high_resolution_clock::now();
    preprocess(tSeries, means, norms, df, dg);
    tend = chrono::high_resolution_clock::now();
    time_elapsed = tend - tstart;
    cout << "[OK] Preprocess Time:         " << setprecision(numeric_limits<double>::digits10 + 2) <<
            time_elapsed.count() << " seconds." << endl;
    /***********************************************/

    /******************** SCAMP ********************/
    cout << "[>>] Performing SCAMP..." << endl;
    tstart = chrono::high_resolution_clock::now();
    scamp(tSeries, means, norms, df, dg, profile, profileIndex);
    tend = chrono::high_resolution_clock::now();
    time_elapsed = tend - tstart;
    cout << "[OK] SCAMP Time:              " << setprecision(numeric_limits<double>::digits10 + 2) <<
            time_elapsed.count() << " seconds." << endl;
    /***********************************************/
 
   
    cout << "[>>] Saving result: " << outfilename << " ..." << endl;
    fstream resultFile(outfilename, ios_base::out);

    for (ITYPE i = 0; i < profileLength; i++) {
       resultFile << setprecision(numeric_limits<double>::digits10 + 2) << (double) sqrt(2 * windowSize * (1 - ((double)profile[i]))) <<
              " " << profileIndex[i] << endl;
    }
    resultFile.close();

    cout << endl;
  } catch (exception &e) {
    cout << "Exception: " << e.what() << endl;
  }
}
