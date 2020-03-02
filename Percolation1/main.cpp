#include<iostream>
#include"mException.h"
#include"fracture2D.h"
#include"mFstream.h"
#include"global.h"
//#include"/Users/Hossein/source/repos/Percolation1/Linear_Solver_GPU/Linear_Solver_GPU.h"
//#include"GPU_SOLVER.h"
using namespace std;

//#include <Eigen/Sparse>


//using Eigen::SparseMatrix;



int main() {

	char run;
	cout << "Hello My Program!!!" << endl;
	double p = 0.7;
	//double lmin = 1.0;
	double lf = 1.0;
	double lmax = 10.0;
	double a = 3; 
	double PR = 5.0;
	double PL = 15.0;
	double tethaMin = -90.0 / 180.0*MY_pi;
	double tethaMax = 90.0 / 180.0*MY_pi; 
	double Lx = 10;
	double Ly = 10;
	int fileNum = 0;
	bool connected;
	string summary_address = "D:/Data/summary/summary.txt";
	result sam = { 0 };
	mFstream summary(sam, summary_address);
	
	for (int i = 0; i < 1000; i++) {
		std::cout << "i = " << i + 1 << endl;
		connected = false;
		int n = 0;
		while (!connected) {
			std::cout << "n = " << n + 1 << endl;
			fileNum++;
			n++;
			string folder = "D:/Data/";
			string fileName = "TEXT_" + to_string(fileNum) + ".txt";
			list<dataStruct> fractureNet;
			fracture2D sample(p, lf, tethaMin, tethaMax, Lx, Ly, folder, fileName, PR, PL);
			//fracture2D sample(p, lmin, lmax, a, tethaMin, tethaMax, Lx, Ly);
			connected = sample.getConnected();
			if (connected) {
				sample.RunPercolation();
				result results = sample.getResults();
				results.iteration = n;
				summary.printResults(results, summary_address, fileNum, i+1);
			}

			//sample.~fracture2D();
		}

	}

	//sample.get_initialData(fractureNet);
	//read_write.writeInitaialState(fractureNet);
	

	//GPU_SOLVER GPU_obj;
	


	std::cout << "End Run!" << std::endl;
	cin >> run;
}