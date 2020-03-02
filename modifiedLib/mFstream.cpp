#include "mFstream.h"



mFstream::mFstream()
{
	ofile.open("initial state", fstream::out);
	clusterText.open("cluster", fstream::out);
	//Exceptions.open("Exceptions.txt", fstream::out);
}

mFstream::mFstream(result results, string address) {

	summary.open(address, ios::out);
	summary << "N" << setw(15) << "Lx" << setw(12) << "Ly" << setw(20) << "Frac_Num" << setw(15) << "occupancy_p" << setw(17)
		<< "percolation_P" << setw(15) << "B" << setw(20) << "D" << setw(15) << "Bf" << setw(15) << "Df" << setw(15) << "K1" <<
		setw(16) << "Kb" << setw(20) << "cluster1" << setw(15) << "cluster2" << setw(15) << "intersections" << setw(15)
		<< "iteration" << setw(10) << "file_Num" << setw(11) << "Pc" << setw(17) << "Q_b" << setw(17) << "Q_b_f" << endl;
}


mFstream::~mFstream()
{
}

void mFstream::writeInitaialState(list<dataStruct> initialState) {

	list<dataStruct>::iterator beginning = initialState.begin();
	register list<dataStruct>::iterator it;
	register int i;



	ofile << "N" << setw(18) << "X center" << setw(25) << "Y center" << setw(25) << "Length" << setw(33) << "Orientarion (Raidan)" << endl;
	i = 0;
	
	for (it = beginning; it != initialState.end(); it++) {
		++i;

		ofile << i << setprecision(10) << setw(20 - int(ceil(log10(i) + 1e-6))) << it->xCenter << setw(25) << it->yCenter << setw(25) << to_string(it->length) << setw(25) << to_string(it->orientation) << endl;
	}
		
}

void mFstream::write_p_c(connectionSPMAT it, int N) {
	
	if (N == 1)
	ofile << "N" << setw(18) << "X center" << setw(25) << "Y center" << endl;

	ofile << N << setprecision(10) << setw(20 - int(ceil(log10(N) + 1e-6))) << it.xi << setw(25) << it.yi << endl;


}

void mFstream::writeCluster(list<cross> clusterList, list<dataStruct> initialData) {

	ofstream clusterFinal;
	clusterFinal.open("clusterFinal", fstream::out);
	clusterText << "Fracture" << setw(18) << "Cluster Number" << setw(18) << "X1" << setw(18) << "Y1" << setw(18) << "X2" << setw(18) << "Y2" << endl;
	for (auto itCross = clusterList.begin(); itCross != clusterList.end(); itCross++) {
		clusterFinal << itCross->fracNum << setw(20) << itCross->cluster << endl;
	}
	clusterFinal.close();
}

void mFstream::onlineWriteCluster(list<cross> clusterList, int N, list<dataStruct> initialData) {

	ofstream clusterText;
	string fileName = "CLUSTER" + to_string(N);
	clusterText.open(fileName, fstream::out);
	register int i = 0;
	clusterText << "Fracture" << setw(18) << "Cluster Number" << setw(12) << "X1" << setw(18) << "Y1" << setw(18) << "X2" << setw(18) << "Y2" << endl;
	auto itinitialData = initialData.begin();
	for (auto itCross = clusterList.begin(); itCross != clusterList.end(); itCross++) {
		(++i);
		double X1 = itinitialData->xL;
		double Y1 = itinitialData->yL;
		double X2 = itinitialData->xR;
		double Y2 = itinitialData->yR;
		itinitialData++;
		clusterText << itCross->fracNum << setw(20 - int(ceil(log10(i) + 1e-6))) << itCross->cluster << setw(21) << to_string(X1) << setw(18) << to_string(Y1) << setw(18) << to_string(X2) << setw(18) << to_string(Y2) << endl;;
	}

	clusterText.close();
}

int mFstream::FileSearch(char *rSeek) {
	register int i;
	char str[MAX_STRING_LENGTH];

	clear();                 // clear fail and eof bits
	seekg(0, std::ios::beg); // back to the start!
	*str = '\0';
	do {
		i = ReadWord(str);
		if (!strcmp(str, rSeek)) return -1;
	} while (i);
	
	return 0;		//Nothing found
}

int mFstream::ReadWord(char *rWord) {
	char ch;
	register int i = 0;

	*rWord = '\0';
	do {
		get(ch);
		if (eof()) {
			return 0;		//Nothing has been read
		}
	} while ((ch < 33) || (ch > 126));

	while ((ch > 32) && (ch < 127)) {
		*(rWord + i) = ch;
		i++;
		get(ch);
		if (eof()) {
			*(rWord + i) = '\0';
			return 1;		//Read, but end of file also encountered
		}
	}

	*(rWord + i) = '\0';
	return -1;		//Correct execution
}

void mFstream::percolatedOutput(list<cross> &crossSparse, list<dataStruct> &initialData, string address) {
	fstream output;
	output.open(address, ios::out);
	output << "Fracture Number" << setw(15) << "X" << setw(15) << "Y" << endl;
	for (auto itCross = crossSparse.begin(); itCross != crossSparse.end(); itCross++) {
		if (itCross->connected) {
			auto itData = initialData.begin();
			advance(itData, itCross->fracNum - 1);
			output << itCross->fracNum << setw(35) << itData->xL << setw(15) << itData->yL << endl;
			output << itCross->fracNum << setw(35) << itData->xR << setw(15) << itData->yR << endl;
		}
	}
	output.close();
}

void mFstream::networkOutput(list<dataStruct> &initialData, string address) {
	fstream output;

	output.open(address, ios::out);
	output << "FN" << setw(20) << "X" << setw(20) << "Y" << endl;
	int n = 0;
	for (auto itData = initialData.begin(); itData != initialData.end(); itData++) {
		output << ++n << setw(22) << itData->xL << setw(22) << itData->yL << endl;
		output << ++n << setw(22) << itData->xR << setw(22) << itData->yR << endl;
	}
	output.close();
}

void  mFstream::backboneOutput(list<connectionSPMAT> &p_c, string address) {
	fstream output;
	output.open(address, ios::out);
	output << "LN" << setw(20) << "X" << setw(20) << "Y" << endl;
	int n = 0;
	for (auto itpci = p_c.begin(); itpci != p_c.end(); itpci++) {
		if ((itpci->flag) && (!(itpci->deadend))) {
			for (auto itp = itpci->backbone.begin(); itp != itpci->backbone.end(); itp++) {
				list<connectionSPMAT>::iterator itpcj = p_c.begin();
				advance(itpcj, *itp - 1);
				output << itpci->pointNum << setw(22) << itpci->xi << setw(22) << itpci->yi << endl;
				output << itpci->pointNum << setw(22) << itpcj->xi << setw(22) << itpcj->yi << endl;
			}
		}		
	}
	output.close();
}

void mFstream::ResultOutput(string address) {

}


void  mFstream::backbonefinalOutput(list<connectionSPMAT> &p_c, string address) {
	fstream output;
	output.open(address, ios::out);
	output << "LN" << setw(20) << "X" << setw(20) << "Y" << endl;
	int n = 0;
	for (auto itpci = p_c.begin(); itpci != p_c.end(); itpci++) {
		if (itpci->backboneFinal.size() > 0) {
			for (auto itp = itpci->backboneFinal.begin(); itp != itpci->backboneFinal.end(); itp++) {
				list<connectionSPMAT>::iterator itpcj = p_c.begin();
				advance(itpcj, *itp - 1);
				output << itpci->pointNum << setw(22) << itpci->xi << setw(22) << itpci->yi << endl;
				output << itpci->pointNum << setw(22) << itpcj->xi << setw(22) << itpcj->yi << endl;
			}
		}
	}
	output.close();
}



void mFstream::exceptionsOutput(string msg) {

	Exceptions.open("Exceptions.txt", fstream::out);
	Exceptions << msg << endl;
	Exceptions.close();
}

void mFstream::printResults(result results, string address, int fileNum, int N) {
	summary << N << setw(15) << results.Lx << setw(15) << results.Ly << setw(15) << results.N << setw(15) << results.p << setw(18)
		<< results.P << setw(21) << results.B << setw(18) << results.D << setw(16) << results.Bf << setw(16) << results.Df << setw(15) << results.K <<
		setw(15) << results.Kb << setw(15) << results.cluster1 << setw(15) << results.cluster2 << setw(15) << results.intersections << setw(15)
		<< results.iteration << setw(10) << fileNum << setw(15) << results.Pc << setw(15) << results.Q1 <<setw(18) <<results.Q2 << endl;

}