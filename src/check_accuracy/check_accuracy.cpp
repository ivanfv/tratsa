#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <limits>
#include <chrono>

using namespace std;

#define RANGE 100 //RIC rango para calcular Top-K accuracy

#define K 1000    //RIC top 100

//RIC para hacer la ordenación del perfil necesito una estructura que contenga
//el valor del perfil y los dos índices (uno para la ventana i y otro para la j)
//ya que al ordenar el mprofile se pierde la información de la i

template <class myType>
class profItem {
public:
  myType val;
  int i;
  int j;
};


bool comp_f(profItem<double> a, profItem<double> b) {

  return (a.val > b.val);

}


bool comp_d(profItem<double> a, profItem<double> b) {

  return (a.val > b.val);

}


bool comp_ff(profItem<double> a, profItem<double> b) {

  return (a.val > b.val);

}


void getProfileAccuracy(string outfilename, int profileLength, vector<profItem<double>> &profile_f,
        vector<profItem<double>> &profile_d, vector<profItem<double>> &profile_ff) {
  chrono::high_resolution_clock::time_point tstart, tend;
  chrono::duration<double> telapsed;

  cout << "------------------------------------------------------------" << endl;
  cout << "[>>] Sorting MProfiles ..." << endl;
  tstart = chrono::high_resolution_clock::now();
  sort(profile_f.begin(), profile_f.end(), comp_f);
  sort(profile_d.begin(), profile_d.end(), comp_d);
  sort(profile_ff.begin(), profile_ff.end(), comp_ff);
  tend = chrono::high_resolution_clock::now();
  telapsed = tend - tstart;
  cout << "[OK] Sorting Time: " << setprecision(numeric_limits<double>::digits10 + 2) << telapsed.count() << " seconds." << endl;


	//for(int l = profileLength - 1;  l > profileLength -10; l--) cout << profile_d[l].val << " " <<  profile_d[l].i <<" " <<  profile_d[l].j << endl; 
	
	//cout << endl;

	//for(int l = profileLength - 1;  l > profileLength -10; l--) cout << profile_f[l].val << " " <<  profile_f[l].i <<" " <<  profile_f[l].j << endl; 
	
	//cout << endl;

	//for(int l = 0;  l < 100; l++) cout << profile_ff[l].val << " " <<  profile_ff[l].i <<" " <<  profile_ff[l].j << endl; 


  //RIC calculo el porcentaje de matches del top-K motifs del transpreciso y el float con el double
  int motifs_f = 0, motifs_f_range = 0, motifs_ff = 0, motifs_ff_range = 0;
  int discords_f = 0, discords_f_range = 0, discords_ff = 0, discords_ff_range = 0;

  for (int i = 0; i < K; i++) {
    for (int j = 0; j < K; j++) {
      //RIC me da igual el orden (i,j) (j,i)
      if (((profile_f[i].i == profile_d[j].i) && (profile_f[i].j == profile_d[j].j)) ||
              ((profile_f[i].j == profile_d[j].i) && (profile_f[i].i == profile_d[j].j))) {
        discords_f++;
        break;
      }
    }
    for (int j = 0; j < K; j++) {
      if ((((profile_f[i].i <= profile_d[j].i + RANGE) && (profile_f[i].i >= profile_d[j].i - RANGE)) &&
              ((profile_f[i].j <= profile_d[j].j + RANGE) && (profile_f[i].j >= profile_d[j].j - RANGE))) ||
              (((profile_f[i].j <= profile_d[j].i + RANGE) && (profile_f[i].j >= profile_d[j].i - RANGE)) &&
              ((profile_f[i].i <= profile_d[j].j + RANGE) && (profile_f[i].i >= profile_d[j].j - RANGE)))) {
        discords_f_range++;
        break;
      }
    }
    for (int j = 0; j < K; j++) {
      //RIC me da igual el orden (i,j) (j,i)
      if (((profile_ff[i].i == profile_d[j].i) && (profile_ff[i].j == profile_d[j].j)) ||
              ((profile_ff[i].j == profile_d[j].i) && (profile_ff[i].i == profile_d[j].j))) {
        discords_ff++;
        break;
      }
    }
    for (int j = 0; j < K; j++) {
      if ((((profile_ff[i].i <= profile_d[j].i + RANGE) && (profile_ff[i].i >= profile_d[j].i - RANGE)) &&
              ((profile_ff[i].j <= profile_d[j].j + RANGE) && (profile_ff[i].j >= profile_d[j].j - RANGE))) ||
              (((profile_ff[i].j <= profile_d[j].i + RANGE) && (profile_ff[i].j >= profile_d[j].i - RANGE)) &&
              ((profile_ff[i].i <= profile_d[j].j + RANGE) && (profile_ff[i].i >= profile_d[j].j - RANGE)))) {
        discords_ff_range++;
        break;
      }
    }
}
  for (int i = profileLength - 1; i >= (profileLength - K); i--) {
    for (int j = profileLength - 1; j >= (profileLength - K); j--) {
      if (((profile_f[i].i == profile_d[j].i) && (profile_f[i].j == profile_d[j].j)) ||
              ((profile_f[i].j == profile_d[j].i) && (profile_f[i].i == profile_d[j].j))) {
        motifs_f++;
        break;
      }
    }
    for (int j = profileLength - 1; j >= (profileLength - K); j--) {
      if ((((profile_f[i].i <= profile_d[j].i + RANGE) && (profile_f[i].i >= profile_d[j].i - RANGE)) &&
              ((profile_f[i].j <= profile_d[j].j + RANGE) && (profile_f[i].j >= profile_d[j].j - RANGE))) ||
              (((profile_f[i].j <= profile_d[j].i + RANGE) && (profile_f[i].j >= profile_d[j].i - RANGE)) &&
              ((profile_f[i].i <= profile_d[j].j + RANGE) && (profile_f[i].i >= profile_d[j].j - RANGE)))) {
        motifs_f_range++;
        break;
      }
    }
    for (int j = profileLength - 1; j >= (profileLength - K); j--) {
      if (((profile_ff[i].i == profile_d[j].i) && (profile_ff[i].j == profile_d[j].j)) ||
              ((profile_ff[i].j == profile_d[j].i) && (profile_ff[i].i == profile_d[j].j))) {
        motifs_ff++;
        break;
      }
    }
    for (int j = profileLength - 1; j >= (profileLength - K); j--) {
      if ((((profile_ff[i].i <= profile_d[j].i + RANGE) && (profile_ff[i].i >= profile_d[j].i - RANGE)) &&
              ((profile_ff[i].j <= profile_d[j].j + RANGE) && (profile_ff[i].j >= profile_d[j].j - RANGE))) ||
              (((profile_ff[i].j <= profile_d[j].i + RANGE) && (profile_ff[i].j >= profile_d[j].i - RANGE)) &&
              ((profile_ff[i].i <= profile_d[j].j + RANGE) && (profile_ff[i].i >= profile_d[j].j - RANGE)))) {
        motifs_ff_range++;
        break;
      }
    }
}

  
  cout << "[>>] Top-" << K << " Accuracy: " << endl;
  cout << "Motifs (float): \t" << setprecision(3) << (motifs_f * 100.0) / K << "%" << endl;
  cout << "Motifs (float)(\u00B1" << RANGE << "): \t" << setprecision(3) << (motifs_f_range * 100.0) / K << "%" << endl;
  cout << "Discs. (float): \t" << setprecision(3) << (discords_f * 100.0) / K << "%" << endl;
  cout << "Discs. (float)(\u00B1" << RANGE << "): \t" << setprecision(3) << (discords_f_range * 100.0) / K << "%" << endl;

  cout << "Motifs (tran.):    \t" << setprecision(3) << (motifs_ff * 100.0) / K << "%" << endl;
  cout << "Motifs (tran.)(\u00B1" << RANGE << "): \t" << (motifs_ff_range * 100.0) / K << "%" << endl;
  cout << "Discs. (tran.): \t" << (discords_ff * 100.0) / K << "%" << endl;
  cout << "Discs. (tran.)(\u00B1" << RANGE << "): \t" << (discords_ff_range * 100.0) / K << "%" << endl;

  /*--------------------------------------------------------------------------*/
  // RIC imprimo al final del fichero
  fstream statsFile(outfilename, ios_base::app);
  //RIC imprimo el tiempo en primer lugar
  statsFile << "#Qsort Time (s): float+double+transp" << endl;
  statsFile << telapsed.count() << endl;
  statsFile << "#Top-" << K << " discords float (val,i,j)" << endl;
  int i;
  for (i = 0; i < K - 1; i++) statsFile << profile_f[i].val << ",";
  statsFile << profile_f[i].val << endl;
  for (i = 0; i < K - 1; i++) statsFile << profile_f[i].i << ",";
  statsFile << profile_f[i].i << endl;
  for (i = 0; i < K - 1; i++) statsFile << profile_f[i].j << ",";
  statsFile << profile_f[i].j << endl;
  statsFile << "#Top-" << K << " discords double (val,i,j)" << endl;
  for (i = 0; i < K - 1; i++) statsFile << profile_d[i].val << ",";
  statsFile << profile_d[i].val << endl;
  for (i = 0; i < K - 1; i++) statsFile << profile_d[i].i << ",";
  statsFile << profile_d[i].i << endl;
  for (i = 0; i < K - 1; i++) statsFile << profile_d[i].j << ",";
  statsFile << profile_d[i].j << endl;
  statsFile << "#Top-" << K << " discords transp (val,i,j)" << endl;
  for (i = 0; i < K - 1; i++) statsFile << profile_ff[i].val << ",";
  statsFile << profile_ff[i].val << endl;
  for (i = 0; i < K - 1; i++) statsFile << profile_ff[i].i << ",";
  statsFile << profile_ff[i].i << endl;
  for (i = 0; i < K - 1; i++) statsFile << profile_ff[i].j << ",";
  statsFile << profile_ff[i].j << endl;
  //RIC discords
  statsFile << "#Top-" << K << " motifs float (val,i,j)" << endl;
  for (i = profileLength - 1; i >= (profileLength - K + 1); i--) statsFile << profile_f[i].val << ",";
  statsFile << profile_f[i].val << endl;
  for (i = profileLength - 1; i >= (profileLength - K + 1); i--) statsFile << profile_f[i].i << ",";
  statsFile << profile_f[i].i << endl;
  for (i = profileLength - 1; i >= (profileLength - K + 1); i--) statsFile << profile_f[i].j << ",";
  statsFile << profile_f[i].j << endl;
  statsFile << "#Top-" << K << " motifs double (val,i,j)" << endl;
  for (i = profileLength - 1; i >= (profileLength - K + 1); i--) statsFile <<  profile_d[i].val << ",";
  statsFile << profile_d[i].val << endl;
  for (i = profileLength - 1; i >= (profileLength - K + 1); i--) statsFile << profile_d[i].i << ",";
  statsFile << profile_d[i].i << endl;
  for (i = profileLength - 1; i >= (profileLength - K + 1); i--) statsFile << profile_d[i].j << ",";
  statsFile << profile_d[i].j << endl;
  statsFile << "#Top-" << K << " motifs flexfloat (val,i,j)" << endl;
  for (i = profileLength - 1; i >= (profileLength - K + 1); i--) statsFile << profile_ff[i].val << ",";
  statsFile << profile_ff[i].val << endl;
  for (i = profileLength - 1; i >= (profileLength - K + 1); i--) statsFile << profile_ff[i].i << ",";
  statsFile << profile_ff[i].i << endl;
  for (i = profileLength - 1; i >= (profileLength - K + 1); i--) statsFile << profile_ff[i].j << ",";
  statsFile << profile_ff[i].j << endl;
  statsFile << "#Top-" << K << " accuracy w.r.t double: Motifs float,float\u00B1" << RANGE << ",transp,transp\u00B1" << RANGE << endl;
  statsFile << setprecision(3) << (motifs_f * 100.0) / K << "," << setprecision(3) << (motifs_f_range * 100.0) / K << "," <<
          setprecision(3) << (motifs_ff * 100.0) / K << "," << (motifs_ff_range * 100.0) / K << endl;
  statsFile << "#Top-" << K << " accuracy w.r.t double: Discords float,float\u00B1" << RANGE << ",tranps,transp\u00B1" << RANGE << endl;
  statsFile << setprecision(3) << (discords_f * 100.0) / K << "," << (discords_f_range * 100.0) / K << "," <<
          setprecision(3) << (discords_ff * 100.0) / K << "," << (discords_ff_range * 100.0) / K << endl;
  statsFile.close();
}





int main(int argc, char* argv[]) {

    unsigned profileLength;
	string doublefilename = argv[1];
	string floatfilename  = argv[2];
	string transpfilename = argv[3];
	string outfilename    = argv[4];

    // Creation of time meassure structures
    chrono::high_resolution_clock::time_point tstart, tend;
   	chrono::duration<double> time_elapsed, mp_time_f, mp_time_d, mp_time_ff;

	vector<profItem<double>> mp_double;
	vector<profItem<double>>  mp_float;
	vector<profItem<double>>  mp_transp;

	/* ------------------------------------------------------------------ */
	/* Read time series file */
	tstart = chrono::high_resolution_clock::now();

	fstream doubleFile(doublefilename, ios_base::in);
	fstream floatFile(floatfilename, ios_base::in);
	fstream transpFile(transpfilename, ios_base::in);
	
	profileLength = 0;

	double value;
	unsigned index;

	profItem<double> tmp_val_double;
	profItem<double>  tmp_val_float;
	profItem<double>  tmp_val_transp;

 	std::string line;
 	while (std::getline(doubleFile, line)) {
		std::stringstream ss(line);

		ss >> value;
		ss >> index;

		tmp_val_double.val = (double) value;
		tmp_val_double.i   = index;
		tmp_val_double.j   = profileLength;

	      	mp_double.push_back(tmp_val_double);
	      	profileLength++;
	    }
	    	    doubleFile.close();

	    profileLength = 0;

 	while (std::getline(floatFile, line)) {
	
		std::stringstream ss(line);
	    
		ss >> value;
		ss >> index;
	
		tmp_val_float.val = (double) value;
		tmp_val_float.i   = index;
		tmp_val_float.j   = profileLength;
	      	mp_float.push_back(tmp_val_float);

	      	profileLength++;
	    }

	    floatFile.close();
	    profileLength = 0;

 	while (std::getline(transpFile, line)) {
	
		std::stringstream ss(line);
	    
		ss >> value;
		ss >> index;
	
		tmp_val_transp.val = (double) value;
		tmp_val_transp.i   = index;
		tmp_val_transp.j   = profileLength;
	      	mp_transp.push_back(tmp_val_transp);

	      	profileLength++;
	    }

	    transpFile.close();

	cout << "[>>] Profile Length: " << profileLength << endl;


	    tend = chrono::high_resolution_clock::now();
	    time_elapsed = tend - tstart;
	    cout << "[OK] Read File Time: " << setprecision(numeric_limits<double>::digits10 + 2) << time_elapsed.count() << " seconds." << endl;


	getProfileAccuracy(outfilename, profileLength, mp_float, mp_double, mp_transp);

	return 0;
}


