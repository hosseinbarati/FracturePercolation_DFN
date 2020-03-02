#include "fracture2D.h"



fracture2D::fracture2D(int N) {



}

fracture2D::fracture2D()
{
}

fracture2D::fracture2D(double p, double lmin, double lmax, double a, double tethaMin, double tethaMax, double Lx, double Ly, string folder, string fileName, double Pr, double Pl)
{
	numFrBac = 0;
	PR = Pr;
	PL = Pl;
	nnz = 0;
	string address = folder + fileName;
	mFstream readFileR;
	//try
	//{
	char str[MAX_STRING_LENGTH], str1[MAX_STRING_LENGTH];
	readFileR.open(address, ios::in);
	if (readFileR.is_open()) {
		cout << "reading the file: " << endl << address << endl << "...................." << endl;
	}
	else
	{
		cout << "Error while reading: /n" << endl;
		cout << address;
	}
	char Xc_s[] = "Xc";
	if (!readFileR.FileSearch(Xc_s)) mException("Number of fractures Not Found");
	char Yc_s[] = "Xc";
	if (!readFileR.FileSearch(Yc_s)) mException("Number of fractures Not Found");
	char Orientation_s[] = "Orientation";
	if (!readFileR.FileSearch(Orientation_s)) mException("Number of fractures Not Found");


	//}
	//catch (mException th)
	//{

	//}
	

	totalPointNum = 0;
//	register int i;
	this->p = p;
	this->lmin = lmin;
	this->a = a;
	this->tethaMin = tethaMin;
	this->tethaMax = tethaMax;
	this->Lx = Lx;
	this->Ly = Ly;
	aveProjlX = 0;
	a_ex = 2 / MY_pi;
	this->lmax = lmax;

	list<dataStruct>::iterator lastElement;
	list<dataStruct>::iterator it;
	list<dataStruct>::iterator ending;
	list<connectionSPMAT>::iterator beginningP_c = p_c.begin(); // size = 2*N
	list<connectionSPMAT>::iterator itP_c;
	N = 0;
	mFstream writeOBJ;
	NCluster = 0;
	connected = false;
	clusterConnection tempClCon;
	do {
		NCluster++;
		dataStruct tempInitialData;
		connectionSPMAT tempP_cL, tempP_cR;
		++N;
		if (!readFileR.ReadWord(str)) mException("Incorrect value for Xc in the data file.");
		tempInitialData.xCenter = atof(str);

		if (!readFileR.ReadWord(str)) mException("Incorrect value for Yc in the data file.");
		tempInitialData.yCenter = atof(str);

		if (!readFileR.ReadWord(str)) mException("Incorrect value for Orientation in the data file.");
		tempInitialData.orientation = atof(str);

		cross tempCross;
		tempCross.cluster = NCluster;
		tempCross.fracNum = N;
		tempCross.connected = false;
		tempClCon.LR.first = false;
		tempClCon.LR.second = false;
		cluCon.push_back(tempClCon);
		crossSparse.push_back(tempCross);
		//randomFractureGen(tempInitialData);
		fractureInit(tempInitialData);						// fracture generation
		initialData.push_back(tempInitialData);	/* put fracture props (xCenter, yCenter,orientation,length,xR,yR,xL,yL,slop,c ) in initialData , the fracture number will be counted by N*/

		//duplicate_initialData.push_back(tempInitialData);; // make a copy of initialData
		insP_c(tempP_cL, tempP_cR, tempInitialData);	// this function provides the 2 point of fractures for saving in p_c
		p_c.push_back(tempP_cL);	// p_c holds the connections, size of p_c is 2*N
		p_c.push_back(tempP_cR);	// p_c holds the connections, size of p_c is 2*N

		/*fNXL,...,f1XL,f1XR,fNXR*/
		//writeOBJ.write_p_c(tempP_cL, totalPointNum - 1);
		//writeOBJ.write_p_c(tempP_cR, totalPointNum);

		findingIntersects();


		calOccupancy_p();
		//dataStruct temp = initialData.back();
//		writeOBJ.onlineWriteCluster(crossSparse, N, duplicate_initialData);
		//if (!connected) {
		connected = checkLR(tempInitialData);



		//std::cout << "occupancy probability is : " << occupancy_p << endl;
	} while (occupancy_p < p);
	readFileR.close();
	//std::cout << "number of fracture is:" << N << endl;
	updateFlags();
	updateCross();
}



fracture2D::fracture2D(double p, double lf, double tethaMin, double tethaMax, double Lx, double Ly, string folder, string fileName, double Pr, double Pl)
{
	numFrBac = 0;
	PR = Pr;
	PL = Pl;
	nnz = 0;
	txtFolder = folder;
	txtName = fileName;
	mFstream readFileR;

	char str[MAX_STRING_LENGTH]; // str1[MAX_STRING_LENGTH];
	string address = folder + fileName;
	readFileR.open(address, ios::in);
	if (readFileR.is_open()) {
		cout << "reading the file: " << endl << address << endl << "...................." << endl;
	}
	else
	{
		cout << "Error while reading: /n" << endl;
		cout << address;
	}
	char Xc_s[] = "Xc";
	if (!readFileR.FileSearch(Xc_s)) mException("Number of fractures Not Found");
	char Yc_s[] = "Xc";
	if (!readFileR.FileSearch(Yc_s)) mException("Number of fractures Not Found");
	char Orientation_s[] = "Orientation";
	if (!readFileR.FileSearch(Orientation_s)) mException("Number of fractures Not Found");
	
	
	totalPointNum = 0;
	register int i;
	this->p = p;
	this->lmin = lf;
	this->a = a;
	this->tethaMin = tethaMin;
	this->tethaMax = tethaMax;
	this->Lx = Lx;
	this->Ly = Ly;
	aveProjlX = 0;
	a_ex = 2 / MY_pi;

	list<dataStruct>::iterator lastElement;
	list<dataStruct>::iterator it;
	list<dataStruct>::iterator ending;
	list<connectionSPMAT>::iterator beginningP_c = p_c.begin(); // size = 2*N
	list<connectionSPMAT>::iterator itP_c;
	N = 0;
	mFstream writeOBJ;
	NCluster = 0;
	connected = false;
	clusterConnection tempClCon;
	do {
		NCluster++;
		dataStruct tempInitialData;
		connectionSPMAT tempP_cL, tempP_cR;
		++N;

		if (!readFileR.ReadWord(str)) mException("Incorrect value for Fracture Number in the data file.");
		if (N != atoi(str)) {
			cout << "Line " << N + 1 << endl;
			mException("Error: Check the fracture Number");		
		}

		if (!readFileR.ReadWord(str)) mException("Incorrect value for Xc in the data file.");
		tempInitialData.xCenter = atof(str);

		if (!readFileR.ReadWord(str)) mException("Incorrect value for Yc in the data file.");
		tempInitialData.yCenter = atof(str);

		if (!readFileR.ReadWord(str)) mException("Incorrect value for Orientation in the data file.");
		tempInitialData.orientation = atof(str);

		cross tempCross;
		tempCross.cluster = NCluster;
		tempCross.fracNum = N;
		tempCross.connected = false;	
		tempClCon.LR.first = false;
		tempClCon.LR.second = false;
		cluCon.push_back(tempClCon);
		crossSparse.push_back(tempCross);
		//(tempInitialData);
		fractureInitConstantLength(tempInitialData);						// fracture generation
		initialData.push_back(tempInitialData);	/* put fracture props (xCenter, yCenter,orientation,length,xR,yR,xL,yL,slop,c ) in initialData , the fracture number will be counted by N*/

		//duplicate_initialData.push_back(tempInitialData);; // make a copy of initialData
		insP_c(tempP_cL, tempP_cR, tempInitialData);	// this function provides the 2 point of fractures for saving in p_c
		p_c.push_back(tempP_cL);	// p_c holds the connections, size of p_c is 2*N
		p_c.push_back(tempP_cR);	// p_c holds the connections, size of p_c is 2*N
		
		/*fNXL,...,f1XL,f1XR,fNXR*/
		//writeOBJ.write_p_c(tempP_cL, totalPointNum - 1);
		//writeOBJ.write_p_c(tempP_cR, totalPointNum);
		
		findingIntersects();

		
		calOccupancy_p();
		//dataStruct temp = initialData.back();
//		writeOBJ.onlineWriteCluster(crossSparse, N, duplicate_initialData);
		//if (!connected) {
		if (!connected) {
			connected = checkLR(tempInitialData);
			if (connected) {
				Pc = occupancy_p;
				std::cout << "N = " << N << endl;
				std::cout << "occupancy probability is : " << occupancy_p << endl;
			}
		}
		
	} while (occupancy_p < p);
	std::cout << "number of fracture is:" << N << endl;
	updateFlags();
	updateCross();
}




fracture2D::~fracture2D()
{
	
}


void fracture2D::findingIntersects() {
	list<dataStruct>::iterator beginning = initialData.begin();
	list<dataStruct>::iterator lastEl = beginning;
	advance(lastEl, N - 1);
	register int i;
	i = 0;
	list<dataStruct>::iterator it = beginning;
	double cen_dis;
	double x, y;
	bool add_point = false;
	pair<double, double> fL2(lastEl->xL, lastEl->yL);
	pair<double, double> fR2(lastEl->xR, lastEl->yR);
	for (it = initialData.begin(); it != initialData.end(); it++) {
		i++;
		cen_dis = calLength(it->xCenter, lastEl->xCenter, it->yCenter, lastEl->yCenter);
		if (cen_dis <= ((it->length + lastEl->length) / 2.0)) {
			if ((it->slope == 0) && (lastEl->slope == INFINITY)) {
				x = lastEl->xCenter;
				y = it->yCenter;
			}
			else if ((it->slope == INFINITY) && (lastEl->slope == 0)) {
				x = it->xCenter;
				y = lastEl->yCenter;
			}
			else if ((it->slope == INFINITY) && (lastEl->slope == INFINITY)) {
				continue;
			}
			else if ((it->slope == 0) && (lastEl->slope != 0) && (lastEl->slope != INFINITY)) {
				y = it->yCenter;
				x = (y - lastEl->c) / lastEl->slope;
			}
			else if ((lastEl->slope == 0) && (it->slope != 0) && (it->slope != INFINITY)) {
				y = lastEl->yCenter;
				x = (y - it->c) / it->slope;
			}
			else if ((it->slope == INFINITY) && (lastEl->slope != 0) && (lastEl->slope != INFINITY)) {
				x = it->xCenter;
				y = lastEl->slope * x + lastEl->c;
			}
			else if ((lastEl->slope == INFINITY) && (it->slope != 0) && (it->slope != INFINITY)) {
				x = lastEl->xCenter;
				y = it->slope * x + it->c;
			}
			else if ((lastEl->slope == it->slope) && (it->c != lastEl->c)) {
				continue;
			}
			else if ((lastEl->slope == it->slope) && (it->c == lastEl->c)) {
				if (checkMergeLine(i, N)){
					mergeLine(i, N);
					lastEl = initialData.end();
					beginning = initialData.begin();
				}
				continue;
			}
			else{
				x = (it->c - lastEl->c) / (lastEl->slope - it->slope);
				y = it->slope * x + it->c;
			}
			pair<double, double> fL1(it->xL, it->yL);
			pair<double, double> fR1(it->xR, it->yR);
			pair<double, double> P(x, y);
			add_point = interstionValidation(fL1, fR1, fL2, fR2, P, i, N);
			if (add_point) {
	


				auto tempSize = (int)p_c.size();
				list<connectionSPMAT>::iterator nTemp = p_c.begin();


				list<dataStruct>::iterator n_f = initialData.begin();
				list<dataStruct>::iterator n_fBegin = initialData.begin();
				/////////////////////
				//advance(n_f, i - 1);

				//advance(nTemp, nI - 1);
				n_f = n_fBegin;
				nTemp = p_c.begin();
				advance(n_f, N - 1);
				int nI = n_f->p_n.first;
				advance(nTemp, nI - 1);
				list<subCon>::iterator subCon1 = nTemp->conecction.begin();
				subCon1->pointJ.push_back(totalPointNum);


				nTemp = p_c.begin();
				nI = n_f->p_n.second;
				advance(nTemp, nI - 1);
				subCon1 = nTemp->conecction.begin();
				subCon1->pointJ.push_back(totalPointNum);

				n_f = n_fBegin;
				nTemp = p_c.begin();
				advance(n_f, i - 1);
				nI = n_f->p_n.first;
				advance(nTemp, nI - 1);
				subCon1 = nTemp->conecction.begin();
				subCon1->pointJ.push_back(totalPointNum);

				nTemp = p_c.begin();
				nI = n_f->p_n.second;
				advance(nTemp, nI - 1);
				subCon1 = nTemp->conecction.begin();
				subCon1->pointJ.push_back(totalPointNum);


				//////////////////////////
				{n_f = n_fBegin;
				nTemp = p_c.begin();
				advance(n_f, i - 1);
				nI = n_f->p_n.first;
				advance(nTemp, nI - 1);
				list<subCon>::iterator tempSub = nTemp->conecction.begin();
				list<int>::iterator tempPoint = tempSub->pointJ.begin();
				for (tempPoint = tempSub->pointJ.begin(); tempPoint != tempSub->pointJ.end(); tempPoint++) {
					int findedPoint = *tempPoint;
					if ((findedPoint != n_f->p_n.second) && (findedPoint != totalPointNum)) {
						list<connectionSPMAT>::iterator pointIt = p_c.begin();
						advance(pointIt, findedPoint - 1);
						list<subCon>::iterator inCon = pointIt->conecction.begin();
						for (inCon = pointIt->conecction.begin(); inCon != pointIt->conecction.end(); inCon++) {
							if (inCon->slope == tempSub->slope) {
								inCon->pointJ.push_back(totalPointNum);
								list<connectionSPMAT>::iterator p_c_totalPoint = p_c.begin();
								advance(p_c_totalPoint, totalPointNum - 1);
								list<subCon>::iterator subConTotP = p_c_totalPoint->conecction.begin();
								for (subConTotP = p_c_totalPoint->conecction.begin(); subConTotP != p_c_totalPoint->conecction.end(); subConTotP++) {
									if (inCon->slope == subConTotP->slope) {
										subConTotP->pointJ.push_back(findedPoint);
									}
								}
							}
						}
					}
				}
				}


				{n_f = n_fBegin;
				nTemp = p_c.begin();
				advance(n_f, N - 1);
				nI = n_f->p_n.first;
				advance(nTemp, nI - 1);
				list<subCon>::iterator tempSub = nTemp->conecction.begin();
				list<int>::iterator tempPoint = tempSub->pointJ.begin();
				for (tempPoint = tempSub->pointJ.begin(); tempPoint != tempSub->pointJ.end(); tempPoint++) {
					int findedPoint = *tempPoint;
					if ((findedPoint != n_f->p_n.second) && (findedPoint != totalPointNum)) {
						list<connectionSPMAT>::iterator pointIt = p_c.begin();
						advance(pointIt, findedPoint - 1);
						list<subCon>::iterator inCon = pointIt->conecction.begin();
						for (inCon = pointIt->conecction.begin(); inCon != pointIt->conecction.end(); inCon++) {
							if (inCon->slope == tempSub->slope) {
								inCon->pointJ.push_back(totalPointNum);
								list<connectionSPMAT>::iterator p_c_totalPoint = p_c.begin();
								advance(p_c_totalPoint, totalPointNum - 1);
								list<subCon>::iterator subConTotP = p_c_totalPoint->conecction.begin();
								for (subConTotP = p_c_totalPoint->conecction.begin(); subConTotP != p_c_totalPoint->conecction.end(); subConTotP++) {
									if (inCon->slope == subConTotP->slope) {
										subConTotP->pointJ.push_back(findedPoint);
									}
								}
							}
						}
					}
				}
				}

				list<cross>::iterator endCross = crossSparse.end();
				list<cross>::iterator itCross = crossSparse.begin();
				list<cross>::iterator itCross1 = crossSparse.begin();
				list<cross>::iterator itCross2 = crossSparse.begin();

				advance(itCross1, i - 1);
				itCross1->col.push_back(N);
				
				advance(itCross2, N - 1);
				itCross2->col.push_back(i);
				
				nTemp = p_c.begin();
				advance(nTemp, totalPointNum - 1);
				nTemp->cluster = (itCross1->cluster < itCross2->cluster) ? itCross1->cluster : itCross2->cluster;

				if (itCross1->cluster < itCross2->cluster) {


					int c1 = itCross1->cluster;
					int c2 = itCross2->cluster;
					for (auto itCross4 = itCross; itCross4 != endCross; itCross4++) {

						if (itCross4->cluster == c2) {
							itCross4->cluster = c1;
						}
						if (itCross4->cluster > c2) {
							itCross4->cluster = itCross4->cluster - 1;
						}
					}
					for (auto itPci = p_c.begin(); itPci != p_c.end(); itPci++) {
						if (itPci->cluster == c2) {
							itPci->cluster = c1;
						}
						if (itPci->cluster > c2) {
							itPci->cluster = itPci->cluster - 1;
						}
					}

					list<clusterConnection>::iterator itCluCon = cluCon.begin();
					advance(itCluCon, c2 - 1);
					list<clusterConnection>::iterator itCluCon1 = cluCon.begin();
					advance(itCluCon1, c1 - 1);
					itCluCon1->LR.first = (itCluCon->LR.first) || (itCluCon1->LR.first);
					itCluCon1->LR.second = (itCluCon->LR.second) || (itCluCon1->LR.second);
					cluCon.erase(itCluCon);
					NCluster--;
				}
				else if (itCross1->cluster > itCross2->cluster) {
					
					int c1 = itCross1->cluster;
					int c2 = itCross2->cluster;
					for (auto itCross4 = itCross; itCross4 != endCross; itCross4++) {
						if (itCross4->cluster == c1) {
							itCross4->cluster = c2;
						}
						if (itCross4->cluster > c1) {
							itCross4->cluster = itCross4->cluster - 1;
						}
					}
					for (auto itPci = p_c.begin(); itPci != p_c.end(); itPci++) {
						if (itPci->cluster == c1) {
							itPci->cluster = c2;
						}
						if (itPci->cluster > c1) {
							itPci->cluster = itPci->cluster - 1;
						}
					
					}
					list<clusterConnection>::iterator itCluCon = cluCon.begin();
					advance(itCluCon, c1 - 1);
					list<clusterConnection>::iterator itCluCon2 = cluCon.begin();
					advance(itCluCon2, c2 - 1);
					itCluCon2->LR.first = (itCluCon->LR.first) || (itCluCon2->LR.first);
					itCluCon2->LR.second = (itCluCon->LR.second) || (itCluCon2->LR.second);
					cluCon.erase(itCluCon);
					NCluster--;
				}

				

			}

		}


	}
}
bool fracture2D::interstionValidation(pair<double, double> fL1, pair<double, double> fR1, pair<double, double> fL2, pair<double, double> fR2, pair<double, double> P, int i, int j) {

	list<cross>::iterator it = crossSparse.begin();
	list<cross>::iterator beginning = crossSparse.begin();
	list<connectionSPMAT>::iterator itP_c1 = p_c.begin(); // size = 2*N
	list<connectionSPMAT>::iterator itP_c2 = p_c.begin();
	list<connectionSPMAT>::iterator beginning_itPc = p_c.begin();
	//list<connectionSPMAT>::iterator itP_cEnd = --p_c.end();
	list<dataStruct>::iterator itDuIniData1 = initialData.begin();
	list<dataStruct>::iterator itDuIniData2 = initialData.begin();

	subCon tempSubcon;
	connectionSPMAT tempP_c;
	dataStruct tempDatastruct;
//	int fNum;
	if ((P.first >= fL1.first) && (P.first <= fR1.first))
		if ((P.first >= fL2.first) && (P.first <= fR2.first))
			if ((P.second >= min(fL1.second, fR1.second)) && (P.second <= max(fL1.second, fR1.second)))
				if ((P.second >= min(fL2.second, fR2.second)) && (P.second <= max(fL2.second, fR2.second))) {

					++totalPointNum;
					tempP_c.pointNum = totalPointNum;
					tempP_c.cluster_fb = -1;


					
					itP_c2 = p_c.begin();


					//advance(itP_c, N - 1);
					// add new point
					tempP_c.xi = P.first;
					tempP_c.yi = P.second;
					// add new point connection N+1
					//////////////////////////////////////
					{subCon temp2;
					itDuIniData2 = itDuIniData1;
					advance(itDuIniData2, N - 1);
					temp2.pointJ.push_back(itDuIniData2->p_n.first);
					temp2.pointJ.push_back(itDuIniData2->p_n.second);	 //add second fracture connection to the new point
					temp2.c = itDuIniData2->c;
					temp2.slope = itDuIniData2->slope;
					temp2.b = itDuIniData2->b;
					tempP_c.conecction.push_back(temp2); }



					/*itP_c1 = itP_c2;
					advance(itP_c1, itDuIniData2->p_n.first - 1);*/
					/*auto it = itP_c1->conecction.begin();
					it->pointJ.push_back(totalPointNum);*/
					
					/*for (auto it = itP_c1->conecction.begin(); it != itP_c1->conecction.end(); it++) {
						if (it->slope == temp2.slope) {
							it->pointJ.push_back(totalPointNum);
						}
					}
					itP_c1 = itP_c2;
					advance(itP_c1, itDuIniData2->p_n.second - 1);*/
					/*for (auto it = itP_c1->conecction.begin(); it != itP_c1->conecction.end(); it++) {
						if (it->slope == temp2.slope) {
							it->pointJ.push_back(totalPointNum);
						}
					}
					}*/


					




					////////////////////////////////////

					//itP_c2 = p_c.begin();
					{subCon temp1;
					itDuIniData2 = itDuIniData1;
					advance(itDuIniData2, i - 1);

					temp1.pointJ.push_back(itDuIniData2->p_n.first);
					temp1.pointJ.push_back(itDuIniData2->p_n.second);	 //add first fracture connection to the new point
					temp1.c = itDuIniData2->c;
					temp1.slope = itDuIniData2->slope;
					temp1.b = itDuIniData2->b;
					tempP_c.conecction.push_back(temp1); }

					if (tempP_c.xi == 0) {
						tempP_c.left = true;
					}
					else {
						tempP_c.left = false;
					}

					if (tempP_c.xi == Lx) {
						tempP_c.right = true;
					}
					else {
						tempP_c.right = false;
					}
					tempP_c.deadend = false;
					tempP_c.dangling = false;
					tempP_c.flag = false;
					p_c.push_back(tempP_c);
					/*itP_c1 = itP_c2;
					advance(itP_c1, itDuIniData2->p_n.first - 1);
					for (auto it = itP_c1->conecction.begin(); it != itP_c1->conecction.end(); it++) {
						if (it->slope == temp1.slope) {
							it->pointJ.push_back(totalPointNum);
						}
					}
					itP_c1 = itP_c2;
					advance(itP_c1, itDuIniData2->p_n.second - 1);
					for (auto it = itP_c1->conecction.begin(); it != itP_c1->conecction.end(); it++){
						if (it->slope == temp1.slope) {
							it->pointJ.push_back(totalPointNum);
						}
					}
					}*/
					
					//
					
					// add new join to point connections

					/*list<subCon>::iterator subConIt1 = tempP_c.conecction.begin();
					list<subCon>::iterator subConIt2 = tempP_c.conecction.begin();
					
					for (auto it = subConIt1; it != tempP_c.conecction.end(); it++) {
						auto pointPTR1 = it->pointJ.begin();			 
						for (auto jt = pointPTR1; jt != it->pointJ.end(); jt++) {
							subCon tempSubcon2;
							//tempSubcon2.pointJ.push_back(pointNum);
							//tempSubcon2.c = tempSubcon.
							itP_c2 = itP_c1;
							advance(itP_c2, *jt - 1);
							tempSubcon2.pointJ.push_back(tempP_c.pointNum);
							tempSubcon2.c = it->c;
							tempSubcon2.slope = it->slope;
							itP_c2->conecction.push_back(tempSubcon2);
						}
					}*/
					



					/////


					//
					return true;
				}

	return false;
}

/*void fracture2D::mergingOverlapLines(list<dataStruct>::iterator it, list<dataStruct>::iterator ending) {
	register int i;
	list<dataStruct>::iterator beginning = it;
	for (i = 1; i < N; i++) {
		if ((it->slope == ending->slope) && (it->slope != INFINITY)) {
			if (it->c == ending->c) {

			}
		}


		it++;
	}

}*/

void fracture2D::insP_c(connectionSPMAT &outL, connectionSPMAT &outR, dataStruct &in) {

	subCon temp;
	temp.pointJ.push_back(in.p_n.second);
	temp.c = in.c;
	temp.slope = in.slope;
	temp.b = in.b;
	outL.conecction.push_back(temp);
	outL.xi = in.xL;
	outL.yi = in.yL;
	outL.pointNum = in.p_n.first;
	if (outL.xi == 0) {
		outL.left = true;
	}
	else {
		outL.left = false;
	}
	outL.right = false;
	outL.flag = false;
	outL.dangling = false;
	outL.deadend = false;
	outL.cluster = NCluster;
	outL.cluster_fb = -1;
	//outL.merge = false;
	//outL.merged = false;

	subCon temp2;
	temp2.pointJ.push_back(in.p_n.first);
	temp2.c = in.c;
	temp2.slope = in.slope;
	temp2.b = in.b;
	outR.conecction.push_back(temp2);
	outR.xi = in.xR;
	outR.yi = in.yR;
	outR.pointNum = in.p_n.second;
	if (outR.xi == Lx) {
		outR.right = true;
	}
	else {
		outR.right = false;
	}
	outR.left = false;
	outR.flag = false;
	outR.dangling = false;
	outR.deadend = false;
	//outR.merge = false;
	//outR.merged = false;
	outR.cluster = NCluster;
	outR.cluster_fb = -1;
}

void fracture2D::updateCSPMAT(connectionSPMAT insData) {

	p_c.push_back(insData);
}

void fracture2D::calAveExArea(list<dataStruct>::iterator lastElement) {
	double l = lastElement->length;
	double exArea = calExArea(l);
	a_ex = ((N - 1)*a_ex + exArea) / N;
}

double fracture2D::calExArea(double l) {

	return  l*(2.0/MY_pi);
}

void fracture2D::calOccupancy_p() {

	occupancy_p = 1 - exp(-N / (4 * Lx*Lx)*a_ex);
}

void fracture2D::fractureInit(dataStruct & temp) {
	

	temp.length = lengthPowerLawDist();
	temp.xR = fractureDelX(temp.orientation, temp.length) + temp.xCenter;
	temp.xL = -fractureDelX(temp.orientation, temp.length) + temp.xCenter;
	temp.yR = fractureDelY(temp.orientation, temp.length) + temp.yCenter;
	temp.yL = -fractureDelY(temp.orientation, temp.length) + temp.yCenter;
	temp.p_n.first = ++totalPointNum;
	temp.p_n.second = ++totalPointNum;
	temp.slope = tan(temp.orientation);
	temp.b = 1.0;

	if ((abs(temp.orientation - MY_pi / 2) < 1e-10) || (abs(temp.orientation + MY_pi / 2) < 1e-10)) {
		temp.slope = INFINITY;
	}
	if (temp.slope != INFINITY) {
		temp.c = temp.yR - temp.slope*temp.xR;
	}
	else {
		temp.c = temp.xCenter;
	}
	limitBoundary(temp);

}

void fracture2D::randomFractureGen(dataStruct &temp) {
	temp.orientation = fractureOrientation();
	temp.xCenter = fractureCenterLocation(Lx);
	temp.yCenter = fractureCenterLocation(Ly);
}


void fracture2D::fractureInitConstantLength(dataStruct &temp) {

	temp.length = lmin;
	temp.xR = fractureDelX(temp.orientation, temp.length) + temp.xCenter;
	temp.xL = -fractureDelX(temp.orientation, temp.length) + temp.xCenter;
	temp.yR = fractureDelY(temp.orientation, temp.length) + temp.yCenter;
	temp.yL = -fractureDelY(temp.orientation, temp.length) + temp.yCenter;
	temp.p_n.first = ++totalPointNum;
	temp.p_n.second = ++totalPointNum;
	temp.slope = tan(temp.orientation);
	temp.b = 1e-4;

	//temp.slope = (temp.yR - temp.yL) / (temp.xR - temp.xL);
	if ((abs(temp.orientation - MY_pi / 2) < 1e-10) || (abs(temp.orientation + MY_pi / 2) < 1e-10)) {
		temp.slope = INFINITY;
	}
	if (temp.slope != INFINITY) {
		temp.c = temp.yR - temp.slope*temp.xR;
	}
	else {
		temp.c = temp.xCenter;
	}
	limitBoundary(temp);
}



void fracture2D::limitBoundary(dataStruct &temp) {

	bool change = false;

	if (temp.xR > Lx) {
		temp.xR = Lx;
		temp.yR = temp.slope*temp.xR + temp.c;
		change = true;
	}
	if (temp.xL < 0.0) {
		temp.xL = 0.0;
		temp.yL = temp.slope*temp.xL + temp.c;
		change = true;
	}
	if (temp.yR > Ly) {
		temp.yR = Ly;
		temp.xR = (temp.yR - temp.c) / temp.slope;
		change = true;
	}
	if (temp.yR < 0.0) {
		temp.yR = 0.0;
		temp.xR = (temp.yR - temp.c) / temp.slope;
		change = true;
	}
	if (temp.yL > Ly) {
		temp.yL = Ly;
		temp.xL = (temp.yL - temp.c) / temp.slope;
		change = true;
	}
	if (temp.yL < 0.0) {
		temp.yL = 0.0;
		temp.xL = (temp.yL - temp.c) / temp.slope;
		change = true;
	}


	if ((temp.xL < 0) || (temp.xR > Lx) || (temp.yL < 0) || (temp.yR > Ly)){
		std::cout << "WRONG OPERATIO" << endl;
		int yy;
		std::cin >> yy;
	}
	if (temp.xL < 0.0) {
		temp.xL = 0.0;
		temp.yL = tan(temp.orientation)*temp.xL + temp.c;
		change = true;
	}
	if (temp.yR > Ly) {
		temp.yR = Ly;
		temp.xR = (temp.yR - temp.c) / tan(temp.orientation);
		change = true;
	}

	temp.length = calLength(temp.xL, temp.xR, temp.yL, temp.yR);
	temp.xCenter = (temp.xL + temp.xR) / 2.0;
	temp.yCenter = (temp.yL + temp.yR) / 2.0;

}


double fracture2D::calLength(double x1, double x2, double y1, double y2) {

	return sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
}


double fracture2D::lengthPowerLawDist() {
	
	double prob = uniformDist();
	double l = pow(prob, 1.0 / (-a + 1.0))*lmin;
	return (l > lmax) ? lmax : l;
}

double fracture2D::fractureOrientation() {
	double prob = uniformDist();
	return prob * (tethaMax - tethaMin) - (tethaMax - tethaMin) / 2.0;
}

double fracture2D::fractureCenterLocation(double L) {

	return uniformDist()*L;
}

double fracture2D::fractureDelY(double orientation, double length) {

	return sin(orientation)*length / 2.0;
}

double fracture2D::fractureDelX(double orientation, double length){
	
	return cos(orientation)*length / 2.0;
}

void fracture2D::get_initialData(list<dataStruct> &inputStruc) {
	inputStruc.assign(initialData.begin(), initialData.end());
}

void fracture2D::sysEffLength() {
	
	list<dataStruct>::iterator beginning = initialData.begin();
	for (auto it = beginning; it != initialData.end(); it++) {
		aveProjlX += 2 * it->length / MY_pi;
	}
	effLx = MY_pi * Lx / (2 * aveProjlX);
}

double fracture2D::uniformDist() {
	std::random_device rd;
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	//std::mt19937 generator(rd());
	//std::mt19937_64 generator(rd());
	std::uniform_real_distribution<double> distribution(0.0, 1.0);	
	return distribution(generator);
}



list<connectionSPMAT> fracture2D::get_p_c(list<connectionSPMAT> out_p_c) const {
	out_p_c.assign(p_c.begin(), p_c.end());
	return out_p_c;
}
list<connectionSPMAT> fracture2D::set_p_c(list<connectionSPMAT> in_p_c) {
	p_c.assign(in_p_c.begin(), in_p_c.end());
	return in_p_c;
}

bool fracture2D::checkLR(const dataStruct &tempInitialData) {
	
	bool Lchange = false;
	bool Rchange = false;
	list<cross>::iterator itCross = crossSparse.begin();
	advance(itCross, N - 1);
	if (tempInitialData.xL == 0.0) {
		int c = itCross->cluster;
		list<clusterConnection>::iterator itCl = cluCon.begin();
		advance(itCl, c - 1);
		itCl->LR.first = true;
		Lchange = true;
	}
	if (tempInitialData.xR == Lx) {
		int c = itCross->cluster;
		list<clusterConnection>::iterator itCl = cluCon.begin();
		advance(itCl, c - 1);
		itCl->LR.second = true;
		Rchange = true;
	}
	
	int c = 0;
	for (auto itCl = cluCon.begin(); itCl != cluCon.end(); itCl++) {
		++c;
		if ((itCl->LR.first == true) && (itCl->LR.second == true)) {
			//itCross->connected = true;
			std::cout << "The system is connected at p = " << occupancy_p << endl;
			std::cout << "The connected cluster is" << c << endl;
			
			//int a;
			//std::cin >> a;
			return true;
		}
	}

	return false;
}


/*void fracture2D::mergePoints() {
	std::cout << "START MERGING" << endl;
	for (auto itPi = p_c.begin(); itPi != p_c.end(); itPi++) { 
		if (itPi->flag) {
			auto tempPi = itPi;
			for (auto itPj = ++tempPi; itPj != p_c.end(); itPj++) {
				if (itPj->flag) {
					//p_c_equality(*itPi, *itPj);
					if ((almostEqual(itPi->xi, itPj->xi, 1e-32)) && (almostEqual(itPi->yi, itPj->yi, 1e-32))) {
						std::cout << "EQUAL" << endl;
						mFstream obj;
						obj.exceptionsOutput("Two points are equal, bad luck!");
						mException("Two points are equal, bad luck!");
						exit(0);
						
						////if (!itPj->merge) {
						//	for (auto subConJ = itPj->conecction.begin(); subConJ != itPj->conecction.end(); subConJ++) {
						//		for (auto subConI = itPi->conecction.begin(); subConI = itPi->conecction.end(); subConI++) {
						//			if (subConJ->slope == subConI->slope) {
						//				subConI->pointJ
						//			}
						//			itPi->conecction.push_back(*subConJ);
						//			itPj->merge = true;
						//		}
						//	}
						////}
					}
					//else {
					//	itPj->merge = false;
					//}
				}
			}
		}
	}
	std::cout << "MERGING FINISHED" << endl;

}
*/
bool fracture2D::p_c_equality(connectionSPMAT p_c_1, connectionSPMAT p_c_2) {
	if ((almostEqual(p_c_1.xi, p_c_2.xi, 0.0)) && (almostEqual(p_c_1.yi, p_c_2.yi, 0.0))) {
		return true;
	}
	else {
		return false;
	}
}

/*connectionSPMAT fracture2D::merge_p_c(connectionSPMAT &temp1, connectionSPMAT &temp2) {
	
	for (auto itSubCon2 = temp2.conecction.begin(); itSubCon2 != temp2.conecction.end(); itSubCon2++) {
		temp1.conecction.push_back(*itSubCon2);				// must be checked
	}
	temp2.merged = true;
}*/

void fracture2D::RunPercolation() {
	//networkResult();
	deleteingIsolatedFractures();
	//mergePoints(); //duplicate_initialData
	findNearestPoint();
	percolationProbability();
	//percolationResult();
	//p_c_equality(connectionSPMAT p_c_1, connectionSPMAT p_c_2)
	deletingDeadEndFractures();
	//backboneResult();
	backboneProbability();
	backbone_indexing();
	K1 = calAbsPerm();
	deleting_dangling();
	//deleteingIsolatedFractures_d();
	deletingDeadEndFractures_d();
	finalBackboneClustering();
	bool connectivity = finalConnectedCluster();
	//backbone_finalResult();
	backboneFinalProbability();

	if (connectivity) {
		backbone_indexing_final();
		K2 = calAbsPerm_final();
	}
	else {
		K2 = K1;
	}
	setResults();
}

void fracture2D::deleteingIsolatedFractures() {
	int frac;
	int pf, ps;

	list<dataStruct>::iterator iniData1 = initialData.begin();
	list<dataStruct>::iterator iniData2;
	list<connectionSPMAT>::iterator it_p_c1 = p_c.begin();
	list<connectionSPMAT>::iterator it_p_c2;
	int count = 0;
	for (auto itCross = crossSparse.begin(); itCross != crossSparse.end(); itCross++) {
		int clu = itCross->cluster;;
		auto itCluster = cluCon.begin();
		advance(itCluster, clu - 1);

		if ((itCluster->LR.first)&&(itCluster->LR.second)) {		
			frac = itCross->fracNum;
			iniData2 = iniData1;
			advance(iniData2, frac - 1);
			pf = iniData2->p_n.first;
			it_p_c2 = it_p_c1;
			advance(it_p_c2, pf - 1);
			/*if (!(it_p_c2->merge)) {
				it_p_c2->flag = true;
			}*/
			
			ps = iniData2->p_n.second;
			it_p_c2 = it_p_c1;
			advance(it_p_c2, ps - 1);
			/*if (!(it_p_c2->merge)) {
				it_p_c2->flag = true;
			}*/

		}
	}
	
	

}

/*void fracture2D::deleteingIsolatedFractures_d() {
	int frac;
	int pf, ps;

	list<dataStruct>::iterator iniData1 = initialData.begin();
	list<dataStruct>::iterator iniData2;
	list<connectionSPMAT>::iterator it_p_c1 = p_c.begin();
	list<connectionSPMAT>::iterator it_p_c2;
	int count = 0;
	for (auto itCross = crossSparse.begin(); itCross != crossSparse.end(); itCross++) {
		int clu = itCross->cluster;;
		auto itCluster = cluCon.begin();
		advance(itCluster, clu - 1);

		if ((itCluster->LR.first) && (itCluster->LR.second)) {
			frac = itCross->fracNum;
			iniData2 = iniData1;
			advance(iniData2, frac - 1);
			pf = iniData2->p_n.first;
			it_p_c2 = it_p_c1;
			advance(it_p_c2, pf - 1);
			/*if (!(it_p_c2->merge)) {
				it_p_c2->flag = true;
			}

			ps = iniData2->p_n.second;
			it_p_c2 = it_p_c1;
			advance(it_p_c2, ps - 1);
			if (!(it_p_c2->merge)) {
				it_p_c2->flag = true;
			}
	
		}
	}



}*/



void fracture2D::findNearestPoint() {

	int n_of_p;
	double disOld;
	double disNew;
	int n_m;
	int n_p;
	for (auto itPci = p_c.begin(); itPci != p_c.end(); itPci++) {
		//auto itNearP = itPci->nearestPoint.begin();
		n_of_p = 0;
		if (itPci->flag) {
			for (auto itSubCon = itPci->conecction.begin(); itSubCon != itPci->conecction.end(); itSubCon++) {
				disOld = 1e+10;
				disNew = 1e+10;
				n_m = 0;
				n_p = 0;
				for (auto itPj = itSubCon->pointJ.begin(); itPj != itSubCon->pointJ.end(); itPj++) {
					auto itPcj = p_c.begin();
					advance(itPcj, *itPj - 1);
					if ((itPci->xi < itPcj->xi)) {
						disNew = (itPci->xi - itPcj->xi)*(itPci->xi - itPcj->xi) + (itPci->yi - itPcj->yi)*(itPci->yi - itPcj->yi);
						if ((disNew < disOld) && (n_m == 0)) {
							n_m = 1;
							disOld = disNew;
							itPci->nearestPoint.push_back(itPcj->pointNum);
							itPci->backbone.push_back(itPcj->pointNum);
							itPci->backboneFinal.push_back(itPcj->pointNum);

							itPci->nearestPoint_b.push_back(itSubCon->b);
							itPci->backbone_b.push_back(itSubCon->b);
							itPci->backboneFinal_b.push_back(itSubCon->b);
						}
						else if ((disNew < disOld) && (itPci->xi < itPcj->xi) && (n_m == 1))
						{
							disOld = disNew;
							itPci->nearestPoint.pop_back();
							itPci->backbone.pop_back();
							itPci->backboneFinal.pop_back();
							itPci->nearestPoint.push_back(itPcj->pointNum);
							itPci->backbone.push_back(itPcj->pointNum);
							itPci->backboneFinal.push_back(itPcj->pointNum);
						}
					}
				}
				disOld = 1e+10;
				disNew = 1e+10;
				n_m = 0;
				n_p = 0;
				for (auto itPj = itSubCon->pointJ.begin(); itPj != itSubCon->pointJ.end(); itPj++) {
					auto itPcj = p_c.begin();
					advance(itPcj, *itPj - 1);
					if ((itPci->xi > itPcj->xi)) {
						disNew = (itPci->xi - itPcj->xi)*(itPci->xi - itPcj->xi) + (itPci->yi - itPcj->yi)*(itPci->yi - itPcj->yi);
						if ((disNew < disOld) && (n_m == 0)) {
							n_m = 1;
							disOld = disNew;
							itPci->nearestPoint.push_back(itPcj->pointNum);
							itPci->backbone.push_back(itPcj->pointNum);
							itPci->backboneFinal.push_back(itPcj->pointNum);

							itPci->nearestPoint_b.push_back(itSubCon->b);
							itPci->backbone_b.push_back(itSubCon->b);
							itPci->backboneFinal_b.push_back(itSubCon->b);
						}
						else if ((disNew < disOld) && (itPci->xi > itPcj->xi) && (n_m == 1))
						{
							disOld = disNew;
							itPci->nearestPoint.pop_back();
							itPci->backbone.pop_back();
							itPci->backboneFinal.pop_back();
							itPci->nearestPoint.push_back(itPcj->pointNum);
							itPci->backbone.push_back(itPcj->pointNum);
							itPci->backboneFinal.push_back(itPcj->pointNum);
						}
					}

				}
			}
		}
		

	}

}

void fracture2D::percolationProbability() {

	parameter_P = 0;
	S_OP_l = 0;
	S_PerP_l = 0;

	for (auto itPci = p_c.begin(); itPci != p_c.end(); itPci++) {
		if (itPci->flag) {
			for (auto pointsList = itPci->nearestPoint.begin(); pointsList != itPci->nearestPoint.end(); pointsList++) {
				if (*pointsList > itPci->pointNum) {
					auto itPcj = p_c.begin();
					advance(itPcj, *pointsList - 1);
					S_PerP_l += (itPci->xi - itPcj->xi)*(itPci->xi - itPcj->xi) + (itPci->yi - itPcj->yi)*(itPci->yi - itPcj->yi);
				}
			}
		}
	}


	for (auto itFrac = initialData.begin(); itFrac != initialData.end(); itFrac++) {
		S_OP_l += (itFrac->xL - itFrac->xR)*(itFrac->xL - itFrac->xR) + (itFrac->yL - itFrac->yR)*(itFrac->yL - itFrac->yR);
	}

	S_PerP_l /= 2;
	parameter_P = S_PerP_l / S_OP_l;
	//std::cout << "Percolation Probability" << parameter_P << std::endl;

}

void fracture2D::deletingDeadEndFractures() {

	for (auto itPci = p_c.begin(); itPci != p_c.end(); itPci++) {
			findingDeadends(itPci);
	}
		/*if (itPci->flag) {
			if (int(itPci->backbone.size()) < 2) {
				itPci->deadend = true;
				for (auto itPcj = p_c.begin(); itPcj != p_c.end(); itPcj++) {
					if (itPcj->flag) {
						for (auto nearP = itPcj->backbone.begin(); nearP != itPcj->backbone.end(); nearP++) {
							if (itPci->flag) {
								if (*nearP == itPci->pointNum) {
									cout << "Point " << *nearP << " deleted" << endl;
									itPcj->deadEndP.push_back(*nearP);
									itPcj->backbone.erase(nearP);
									break;
								}
							}
						}
					}
				}
			}
		}
	}*/

}


void fracture2D::deletingDeadEndFractures_d() {

	for (auto itPci = p_c.begin(); itPci != p_c.end(); itPci++) {
		findingDeadends_d(itPci);
	}
	/*if (itPci->flag) {
		if (int(itPci->backbone.size()) < 2) {
			itPci->deadend = true;
			for (auto itPcj = p_c.begin(); itPcj != p_c.end(); itPcj++) {
				if (itPcj->flag) {
					for (auto nearP = itPcj->backbone.begin(); nearP != itPcj->backbone.end(); nearP++) {
						if (itPci->flag) {
							if (*nearP == itPci->pointNum) {
								cout << "Point " << *nearP << " deleted" << endl;
								itPcj->deadEndP.push_back(*nearP);
								itPcj->backbone.erase(nearP);
								break;
							}
						}
					}
				}
			}
		}
	}
}*/

}





void fracture2D::findingDeadends(list<connectionSPMAT>::iterator &tempPci) {
	if ((tempPci->flag) && (!(tempPci->deadend))) {
		if ((int(tempPci->backbone.size()) == 1) && (!(tempPci->left)) && (!(tempPci->right))) {
			int Npoint = tempPci->backbone.back();
			double N_b = tempPci->backbone_b.back();

			tempPci->backbone.pop_back();
			tempPci->backboneFinal.pop_back();
			tempPci->backbone_b.pop_back();
			tempPci->backboneFinal_b.pop_back();

			tempPci->deadEndP.push_back(Npoint);
			tempPci->deadEndP_b.push_back(N_b);

			auto itPci = p_c.begin();
			advance(itPci, Npoint - 1);
			auto itp1 = find(itPci->backbone.begin(), itPci->backbone.end(), tempPci->pointNum);
			itPci->backbone.erase(itp1);
			auto it_b1 = find(itPci->backbone_b.begin(), itPci->backbone_b.end(), N_b);
			itPci->backbone_b.erase(it_b1);
			auto itp2 = find(itPci->backboneFinal.begin(), itPci->backboneFinal.end(), tempPci->pointNum);
			itPci->backboneFinal.erase(itp2);
			auto it_b2 = find(itPci->backboneFinal_b.begin(), itPci->backboneFinal_b.end(), N_b);
			itPci->backboneFinal_b.erase(it_b2);



			itPci->deadEndP.push_back(tempPci->pointNum);
			itPci->deadEndP_b.push_back(N_b);
			tempPci->deadend = true;
			//std::cout << "point " << tempPci->pointNum << " deleted from " << itPci->pointNum << endl;
			findingDeadends(itPci);
		}
		else if ((int(tempPci->backbone.size()) == 0) && (!tempPci->deadend)) {
			tempPci->deadend = true;
		}
	}
}


void fracture2D::findingDeadends_d(list<connectionSPMAT>::iterator &tempPci) {
	register int n;
	if ((tempPci->flag) && (!(tempPci->deadend)) && (!tempPci->dangling)) {
		if ((int(tempPci->backboneFinal.size()) == 1) && (!(tempPci->left)) && (!(tempPci->right))) {
			int Npoint = tempPci->backboneFinal.back();
			double N_b = tempPci->backboneFinal_b.back();
			auto itPci = p_c.begin();
			advance(itPci, Npoint - 1);
			

			tempPci->backboneFinal.pop_back();
			tempPci->backboneFinal_b.pop_back();

			tempPci->danglingP.push_back(Npoint);
			tempPci->danglingP_b.push_back(N_b);



			n = 0;
			auto it_b1 = itPci->backboneFinal_b.begin();
			for (auto itp1 = itPci->backboneFinal.begin(); itp1 != itPci->backboneFinal.end(); itp1++) {
				if (*itp1 == tempPci->pointNum) {
					itPci->backboneFinal.erase(itp1);
					advance(it_b1, n);
					itPci->backboneFinal_b.erase(it_b1);
					break;
				}
				n++;
			}

			itPci->danglingP.push_back(tempPci->pointNum);
			itPci->danglingP_b.push_back(N_b);
			//std::cout << "dangling end point " << tempPci->pointNum << " deleted from " << itPci->pointNum << endl;
			findingDeadends_d(itPci);
			tempPci->dangling = true;
		}
		else if ((int(tempPci->backboneFinal.size()) == 0) && (!tempPci->dangling)) {
			tempPci->dangling = true;
		}
	}
}



void fracture2D::updateFlags() {
	//std::cout << "Updating flags...." << std::endl;
	int n = 0;
	for (auto itCluCon = cluCon.begin(); itCluCon != cluCon.end(); itCluCon++) {
		n++;
		if ((itCluCon->LR.first) && (itCluCon->LR.second)) {
			for (auto itPci = p_c.begin(); itPci != p_c.end(); itPci++) {
				if ((itPci->cluster == n)&&(!(itPci->flag))) {
					itPci->flag = true;
				}
			}
		}
	}
	//std::cout << "Flags updating ended." << std::endl;
}

void fracture2D::updateCross() {
	//std::cout << "Updating Cross...." << std::endl;
	int n = 0;
	for (auto itClu = cluCon.begin(); itClu != cluCon.end(); itClu++) {
		n++;
		if ((itClu->LR.first) && (itClu->LR.second)) {
			for (auto itCross = crossSparse.begin(); itCross != crossSparse.end(); itCross++) {
				if ((itCross->cluster == n)) {
					itCross->connected = true;
				}
			}
		}
	}
	//std::cout << "Flags updating ended." << std::endl;
}


void fracture2D::backboneFinalProbability() {

	parameter_B_final = 0;
	parameter_Dang = 0;
	S_PerB_lf = 0;
	S_PerD_lf = 0;
	for (auto itPci = p_c.begin(); itPci != p_c.end(); itPci++) {
		if (itPci->backboneFinal.size() > 0) {
			for (auto pointsList = itPci->backboneFinal.begin(); pointsList != itPci->backboneFinal.end(); pointsList++) {
				if (*pointsList > itPci->pointNum) {
					auto itPcj = p_c.begin();
					advance(itPcj, *pointsList - 1);
					S_PerB_lf += (itPci->xi - itPcj->xi)*(itPci->xi - itPcj->xi) + (itPci->yi - itPcj->yi)*(itPci->yi - itPcj->yi);
				}
			}
		}
	}

	/*for (auto itPci = p_c.begin(); itPci != p_c.end(); itPci++) {
		if (itPci->backboneFinal.size() > 0) {
			for (auto pointsList = itPci->deadEndP.begin(); pointsList != itPci->deadEndP.end(); pointsList++) {
				if (*pointsList > itPci->pointNum) {
					auto itPcj = p_c.begin();
					advance(itPcj, *pointsList - 1);
					S_PerD_lf += (itPci->xi - itPcj->xi)*(itPci->xi - itPcj->xi) + (itPci->yi - itPcj->yi)*(itPci->yi - itPcj->yi);
				}
			}
		}
	}

	for (auto itPci = p_c.begin(); itPci != p_c.end(); itPci++) {
		if (itPci->backboneFinal.size() > 0) {
			for (auto pointsList = itPci->danglingP.begin(); pointsList != itPci->danglingP.end(); pointsList++) {
				if (*pointsList > itPci->pointNum) {
					auto itPcj = p_c.begin();
					advance(itPcj, *pointsList - 1);
					S_PerD_lf += (itPci->xi - itPcj->xi)*(itPci->xi - itPcj->xi) + (itPci->yi - itPcj->yi)*(itPci->yi - itPcj->yi);
				}
			}
		}
	}
*/
	S_PerB_lf /= 2;
	parameter_B_final = S_PerB_lf / S_OP_l;
	parameter_Dang = parameter_P - parameter_B_final;
	//std::cout << "Backbone Probability" << parameter_B << std::endl;
}


void fracture2D::backboneProbability() {

	parameter_B = 0;
	S_PerB_l = 0;
	parameter_D = 0;
	S_PerD_l = 0;
	for (auto itPci = p_c.begin(); itPci != p_c.end(); itPci++) {
		if (itPci->flag) {
			for (auto pointsList = itPci->backbone.begin(); pointsList != itPci->backbone.end(); pointsList++) {
				if (*pointsList > itPci->pointNum) {
					auto itPcj = p_c.begin();
					advance(itPcj, *pointsList - 1);
					S_PerB_l += (itPci->xi - itPcj->xi)*(itPci->xi - itPcj->xi) + (itPci->yi - itPcj->yi)*(itPci->yi - itPcj->yi);
				}
			}
		}
	}

	for (auto itPci = p_c.begin(); itPci != p_c.end(); itPci++) {
		if (itPci->flag) {
			for (auto pointsList = itPci->deadEndP.begin(); pointsList != itPci->deadEndP.end(); pointsList++) {
				if (*pointsList > itPci->pointNum) {
					auto itPcj = p_c.begin();
					advance(itPcj, *pointsList - 1);
					S_PerD_l += (itPci->xi - itPcj->xi)*(itPci->xi - itPcj->xi) + (itPci->yi - itPcj->yi)*(itPci->yi - itPcj->yi);
				}
			}
		}
	}
	S_PerB_l /= 2;
	S_PerD_l /= 2;
	parameter_B = S_PerB_l / S_OP_l;
	parameter_D = S_PerD_l / S_OP_l;
	//std::cout << "Backbone Probability" << parameter_B << std::endl;
}



void fracture2D::setLx(double Lx) {
	this->Lx = Lx;
}
void fracture2D::setLy(double Ly) {
	this->Ly = Ly;
}
void fracture2D::set_lmin(double lmin) {
	this->lmin = lmin;
}
void fracture2D::set_a(double a) {
	this->a = a;
}
void fracture2D::set_initialData(dataStruct temp) {
	this->initialData.push_back(temp);
}
void fracture2D::set_p_c(connectionSPMAT tempP_c) {
	this->p_c.push_back(tempP_c);
}

void fracture2D::setTotalPoint(int n) {
	this->totalPointNum = n;
}

void fracture2D::setNcluster(int n) {
	this->NCluster = n;
}



int fracture2D::getTotalPoint() {
	return totalPointNum;
}int fracture2D::getNcluster() {
	return NCluster;
}


void fracture2D::setConnected(bool val) {
	this->connected = val;
}

bool fracture2D::getConnected() {
	return connected;
}

void fracture2D::setcluCon(clusterConnection &temp) {
	cluCon.push_back(temp);
}

void fracture2D::setCrossSparseCon(cross &temp) {
	crossSparse.push_back(temp);
}

void fracture2D::setp_c(connectionSPMAT &temp) {
	p_c.push_back(temp);
}

void fracture2D::set_a_ex() {
	a_ex = 2 / MY_pi;
}

void fracture2D::setFracNum(int N) {
	this->N = N;
}


void fracture2D::mergeLine(int line1, int line2) {
	list<dataStruct>::iterator dataLine1 = initialData.begin();
	advance(dataLine1, line1 - 1);
	list<dataStruct>::iterator dataLine2 = initialData.begin();
	advance(dataLine2, line2 - 1);

	if (dataLine1->xL > dataLine2->xL) {
		dataLine1->p_n.first = dataLine2->p_n.first;
		dataLine1->xL = dataLine2->xL;
		dataLine1->yL = dataLine2->yL;
	}
	if (dataLine1->xR < dataLine2->xR) {
		dataLine1->p_n.second = dataLine2->p_n.second;
		dataLine1->xR = dataLine2->xR;
		dataLine1->yR = dataLine2->yR;
	}
	dataLine1->xCenter = (dataLine1->xL + dataLine1->xR) / 2;
	dataLine1->yCenter = (dataLine1->yL + dataLine1->yR) / 2;
	dataLine1->length = sqrt((dataLine1->xL - dataLine1->xR)*(dataLine1->xL - dataLine1->xR) + (dataLine1->yL - dataLine1->yR)*(dataLine1->yL - dataLine1->yR));

	initialData.erase(dataLine2);
}
bool fracture2D::checkMergeLine(int line1, int line2) {
	list<dataStruct>::iterator dataLine1 = initialData.begin();
	advance(dataLine1, line1 - 1);
	list<dataStruct>::iterator dataLine2 = initialData.begin();
	advance(dataLine2, line2 - 1);

	if ((dataLine1->xR > dataLine2->xL) && (dataLine1->xR < dataLine2->xR)) {
		return true;
	} else if ((dataLine2->xR > dataLine1->xL) && (dataLine2->xR < dataLine1->xR)) {
		return true;
	}
	
	return false;
}


void fracture2D::percolationResult() {

	mFstream obj;
	string P_address = txtFolder + "Percolation_results/" + txtName;
	obj.percolatedOutput(crossSparse, initialData, P_address);
}

void fracture2D::networkResult() {
	mFstream obj;
	string N_address = txtFolder + "network_results/" + txtName;
	obj.networkOutput(initialData, N_address);
}

void fracture2D::backboneResult() {
	mFstream obj;
	string B_address = txtFolder + "backbone_results/" + txtName;
	obj.backboneOutput(p_c, B_address);
}

void fracture2D::backbone_finalResult() {
	mFstream obj;
	string B_address = txtFolder + "backbone_finalresults/" + txtName;
	obj.backbonefinalOutput(p_c, B_address);
}


void fracture2D::backbone_indexing()
{
	register int i;
	i = 0;
	for (auto itPci = p_c.begin(); itPci != p_c.end(); itPci++) {
		if ((itPci->flag) && (!itPci->deadend) && (!itPci->left) && (!itPci->right)) {
			itPci->i_b = i++;
			nnz += int(itPci->backbone.size());
		}
	}
	numFrBac = i;

}

void fracture2D::backbone_indexing_final()
{
	register int i;
	i = 0;
	nnz_final = 0;
	for (auto itPci = p_c.begin(); itPci != p_c.end(); itPci++) {
		if (itPci->backboneFinal.size()) {
			itPci->i_b_f = i++;
			nnz_final += int(itPci->backboneFinal.size());
		}
	}
	numFrBac_final = i;

}



double fracture2D::calAbsPerm() {
	auto start1 = chrono::high_resolution_clock::now();

	double K = 0;
	SparseMatrix<double> coef(numFrBac, numFrBac);
	//MatrixXd coef(numFrBac, numFrBac);
	int nonZero = (nnz + numFrBac);
	coef.reserve(nonZero);
	VectorXd u(numFrBac);
	VectorXd x(numFrBac);
	/*SparseVectorXd<double> u(numFrBac);
	VectorXd x(numFrBac);*/

	register int n;
	for (n = 0; n < numFrBac; n++) {
		u(n) = 0.0;
		//x(n) = 1.0;
		coef.insert(n, n) = 0.0;
		/*for (z = 0; z < numFrBac; z++) {
			coef(n, z) = 0.0;
		}*/
	}

	n = 0;
	int jp;
	double d = 0;
	double ci = 0;
	double cj = 0;
	double ui = 0;


	for (auto itPci = p_c.begin(); itPci != p_c.end(); itPci++) {
		if ((!itPci->left) && (!itPci->right) && (!itPci->deadend)&&(itPci->flag)) {
			//i = itPci->pointNum;
			//ci = 0;
			n = 0;
			for (auto itBacP = itPci->backbone.begin(); itBacP != itPci->backbone.end(); itBacP++) {
				jp = *itBacP;
				auto itPcj = p_c.begin();
				advance(itPcj, jp - 1);
				d = sqrt((itPci->xi - itPcj->xi)*(itPci->xi - itPcj->xi) + (itPci->yi - itPcj->yi)*(itPci->yi - itPcj->yi));
				auto it_b = itPci->backbone_b.begin();
				advance(it_b, n);
				cj = (*it_b)*(*it_b)*(*it_b) / d;
				//ci += cj;
				n++;
				if ((!itPcj->left) && (!itPcj->right)) {
					//coef.insert(itPci->i_b, itPcj->i_b) = cj;
					coef.insert(itPci->i_b, itPcj->i_b) = -cj;
				}
				else if (itPcj->left) {
					u(itPci->i_b) += cj*PL;
				}
				else if (itPcj->right) {
					u(itPci->i_b) += cj * PR;
				}

				coef.coeffRef(itPci->i_b, itPci->i_b) += cj;
			}
			//coef.insert(itPci->i_b, itPci->i_b) = ci;


		}
	}
	
	auto start2 = chrono::high_resolution_clock::now();
	//Eigen::ConjugateGradient<SparseMatrix<double>> solver;
	//BiCGSTAB<SparseMatrix<double>> solver;
	ConjugateGradient<SparseMatrix<double>> solver;
	solver.compute(coef);
	x = solver.solve(u);
	auto start3 = chrono::high_resolution_clock::now();
	//x = coef.colPivHouseholderQr().solve(u);
	//x = coef.ldlt().solve(u);
	//mFstream file;
	//string address = txtFolder + "u.txt";
	//file.open(address, ios::out);
	//file << u;
	////x = solver.solve(u);
	//mFstream fileAns;
	//string address2 = txtFolder + "ans.txt";
	//fileAns.open(address2, ios::out);
	//fileAns << x;
	//for (n = 0; n < numFrBac; n++)
	//	cout << x(n) << endl;

	double qr = 0.0;
	double ql = 0.0;
	for (auto itPci = p_c.begin(); itPci != p_c.end(); itPci++) {
		if ((itPci->flag) && (!itPci->deadend)) {
			if (itPci->right) {
				n = 0;
				for (auto itBacP = itPci->backbone.begin(); itBacP != itPci->backbone.end(); itBacP++) {
					jp = *itBacP;
					auto itPcj = p_c.begin();
					advance(itPcj, jp - 1);
					d = sqrt((itPci->xi - itPcj->xi)*(itPci->xi - itPcj->xi) + (itPci->yi - itPcj->yi)*(itPci->yi - itPcj->yi));
					auto it_b = itPci->backbone_b.begin();
					advance(it_b, n);
					qr += (*it_b)*(*it_b) / 12 / DarcyToM2 *(*it_b * mToft) / d * (x(itPcj->i_b) - PR);
					n++;
				}
			}
		}
	}

	for (auto itPci = p_c.begin(); itPci != p_c.end(); itPci++) {
		if ((itPci->flag) && (!itPci->deadend)) {
			if (itPci->left) {
				n = 0;
				for (auto itBacP = itPci->backbone.begin(); itBacP != itPci->backbone.end(); itBacP++) {
					jp = *itBacP;
					auto itPcj = p_c.begin();
					advance(itPcj, jp - 1);
					d = sqrt((itPci->xi - itPcj->xi)*(itPci->xi - itPcj->xi) + (itPci->yi - itPcj->yi)*(itPci->yi - itPcj->yi));
					auto it_b = itPci->backbone_b.begin();
					advance(it_b, n);
					ql += (*it_b)*(*it_b) / 12 / DarcyToM2 * (*it_b * mToft) / d * (PL - x(itPcj->i_b));
					n++;
				}
			}
		}
	}

	if (almostEqual(ql, qr, tol_q)) {
		Q_b = ql;
	}
	else {
		mException("ql is not equal to qr!");
	}

	//cout << "Number of equation is " << numFrBac << endl;
	auto start4 = chrono::high_resolution_clock::now();
	auto dt1 = chrono::duration_cast<chrono::microseconds> (start2 - start1).count();
	auto dt2 = chrono::duration_cast<chrono::microseconds> (start3 - start2).count();
	auto dt3 = chrono::duration_cast<chrono::microseconds> (start4 - start3).count();
	//cout << "CREATING SYSTEM of EQUATIONS: " << dt1 << endl;
	//cout << "SOLVING SYSTEM of EQUATIONS: " << dt2 << endl;
	//cout << "VALIDATION: " << dt3 << endl;
	pressureAssignment_backbone(x);

	K = Q_b * Lx / (Ly * (PL - PR));
	return (K/((1e-4)*(1e-4) / 12 / DarcyToM2));
}

double fracture2D::calAbsPerm_final() {
	auto start1 = chrono::high_resolution_clock::now();

	double K = 0;
	//SparseMatrix<double> coef(numFrBac_final, numFrBac_final);
	MatrixXd coef(numFrBac_final, numFrBac_final);
	int nonZero = (nnz_final + numFrBac_final);
	//coef.reserve(nonZero);
	VectorXd u(numFrBac_final);
	VectorXd x(numFrBac_final);
	/*SparseVectorXd<double> u(numFrBac);
	VectorXd x(numFrBac);*/

	register int n, z;
	for (n = 0; n < numFrBac_final; n++) {
		u(n) = 0.0;
		//x(n) = 1.0;
		//coef.insert(n, n) = 0.0;
		for (z = 0; z < numFrBac_final; z++) {
			coef(n, z) = 0.0;
		}
	}

	n = 0;
	int jp;
	double d = 0;
	double ci = 0;
	double cj = 0;
	double ui = 0;


	for (auto itPci = p_c.begin(); itPci != p_c.end(); itPci++) {
		if (itPci->backboneFinal.size() > 0) {
			//i = itPci->pointNum;
			//ci = 0;
			n = 0;
			for (auto itBacP = itPci->backboneFinal.begin(); itBacP != itPci->backboneFinal.end(); itBacP++) {
				jp = *itBacP;
				auto itPcj = p_c.begin();
				advance(itPcj, jp - 1);
				d = sqrt((itPci->xi - itPcj->xi)*(itPci->xi - itPcj->xi) + (itPci->yi - itPcj->yi)*(itPci->yi - itPcj->yi));
				auto it_b = itPci->backboneFinal_b.begin();
				advance(it_b, n);
				cj = (*it_b)*(*it_b)*(*it_b) / d;
				//ci += cj;
				n++;
				if ((!itPcj->left) && (!itPcj->right)) {
					coef(itPci->i_b_f, itPcj->i_b_f) = -cj;
					//coef.insert(itPci->i_b_f, itPcj->i_b_f) = -cj;
				}
				else if (itPcj->left) {
					u(itPci->i_b_f) += cj * PL;
				}
				else if (itPcj->right) {
					u(itPci->i_b_f) += cj * PR;
				}

				coef(itPci->i_b_f, itPci->i_b_f) += cj;
				//coef.coeffRef(itPci->i_b_f, itPci->i_b_f) += cj;
			}
			//coef.insert(itPci->i_b, itPci->i_b) = ci;


		}
	}


	auto start2 = chrono::high_resolution_clock::now();
	//Eigen::ConjugateGradient<SparseMatrix<double>> solver;
	//BiCGSTAB<SparseMatrix<double>> solver;
	//ConjugateGradient<SparseMatrix<double>> solver;
	//solver.compute(coef);
	//x = solver.solve(u);
	auto start3 = chrono::high_resolution_clock::now();
	//x = coef.colPivHouseholderQr().solve(u);
	x = coef.ldlt().solve(u);
	//mFstream file;
	//string address = txtFolder + "u.txt";
	//file.open(address, ios::out);
	//file << u;
	////x = solver.solve(u);
	//mFstream fileAns;
	//string address2 = txtFolder + "ans.txt";
	//fileAns.open(address2, ios::out);
	//fileAns << x;

	double qr = 0.0;
	double ql = 0.0;
	for (auto itPci = p_c.begin(); itPci != p_c.end(); itPci++) {
		if (itPci->backboneFinal.size() > 0) {
			if (itPci->right) {
				n = 0;
				for (auto itBacP = itPci->backboneFinal.begin(); itBacP != itPci->backboneFinal.end(); itBacP++) {
					jp = *itBacP;
					auto itPcj = p_c.begin();
					advance(itPcj, jp - 1);
					d = sqrt((itPci->xi - itPcj->xi)*(itPci->xi - itPcj->xi) + (itPci->yi - itPcj->yi)*(itPci->yi - itPcj->yi));
					auto it_b = itPci->backboneFinal_b.begin();
					advance(it_b, n);
					qr += (*it_b)*(*it_b) / 12 / DarcyToM2 * (*it_b * mToft) / d * (x(itPcj->i_b_f) - PR);
					n++;
				}
			}
		}
	}

	for (auto itPci = p_c.begin(); itPci != p_c.end(); itPci++) {
		if (itPci->backboneFinal.size() > 0) {
			if (itPci->left) {
				n = 0;
				for (auto itBacP = itPci->backboneFinal.begin(); itBacP != itPci->backboneFinal.end(); itBacP++) {
					jp = *itBacP;
					auto itPcj = p_c.begin();
					advance(itPcj, jp - 1);
					d = sqrt((itPci->xi - itPcj->xi)*(itPci->xi - itPcj->xi) + (itPci->yi - itPcj->yi)*(itPci->yi - itPcj->yi));
					auto it_b = itPci->backboneFinal_b.begin();
					advance(it_b, n);
					ql += (*it_b)*(*it_b) / 12 / DarcyToM2 * (*it_b * mToft) / d * (PL - x(itPcj->i_b_f));
					n++;
				}
			}
		}
	}

	if (almostEqual(ql, qr, tol_q)) {
		Q_b_f = ql;
	}
	else {
		mException("ql is not equal to qr!");
	}

	/*cout << "Number of equation is " << numFrBac << endl;
	auto start4 = chrono::high_resolution_clock::now();
	auto dt1 = chrono::duration_cast<chrono::microseconds> (start2 - start1).count();
	auto dt2 = chrono::duration_cast<chrono::microseconds> (start3 - start2).count();
	auto dt3 = chrono::duration_cast<chrono::microseconds> (start4 - start3).count();
	cout << "CREATING SYSTEM of EQUATIONS: " << dt1 << endl;
	cout << "SOLVING SYSTEM of EQUATIONS: " << dt2 << endl;
	cout << "VALIDATION: " << dt3 << endl;*/
	//pressureAssignment_backbone(x);

	K = Q_b_f / (Ly * (PL - PR) / Lx);
	return (K / ((1e-4)*(1e-4) / 12 / DarcyToM2));
}



void fracture2D::pressureAssignment_backbone(VectorXd x) {
	register int n;
	n = 0;
	for (auto itPci = p_c.begin(); itPci != p_c.end(); itPci++) {
		if ((!itPci->left) && (!itPci->right) && (!itPci->deadend) && (itPci->flag)) {
			itPci->P_b = x(n);
			n++;
		}
	}

}

void fracture2D::deleting_dangling() {
	register int n, i;
	n = 0;
	double d;
	double qi;
	for (auto itPci = p_c.begin(); itPci != p_c.end(); itPci++) {
		if ((!itPci->left) && (!itPci->right) && (!itPci->deadend) && (itPci->flag)) {
			n = 0; 
			qi = 0;
			auto itBacP_f = itPci->backboneFinal.begin();
			while (itBacP_f != itPci->backboneFinal.end()) {
				int jp = *itBacP_f;
				auto itPcj = p_c.begin();
				advance(itPcj, jp - 1);
				d = sqrt((itPci->xi - itPcj->xi)*(itPci->xi - itPcj->xi) + (itPci->yi - itPcj->yi)*(itPci->yi - itPcj->yi));
				auto it_b_f = itPci->backboneFinal_b.begin();
				advance(it_b_f, n);
				qi = (*it_b_f)*(*it_b_f) / 12 / DarcyToM2 * (*it_b_f * mToft) / d * (itPci->P_b - itPcj->P_b);
				
				/*if ((qi/Q_b) < tol_q) {
					itPci->danglingP.push_back(jp);
					itPci->danglingP_b.push_back(*it_b_f);
					itPci->backboneFinal.erase(itBacP_f);
					itBacP_f = itPci->backboneFinal.begin();
					advance(itBacP_f, n);
					itPci->backboneFinal_b.erase(it_b_f);
					it_b_f = itPci->backboneFinal_b.begin();
					advance(it_b_f, n);
				}*/

				if (abs((qi / Q_b)) < tol_q) {
					itPci->danglingP.push_back(jp);
					itPci->danglingP_b.push_back(*it_b_f);
					itPci->backboneFinal.erase(itBacP_f);
					itBacP_f = itPci->backboneFinal.begin();
					advance(itBacP_f, n);
					itPci->backboneFinal_b.erase(it_b_f);
				}
				else {
					itBacP_f++;
					n++;
				}

			}
		}
	}

}

void fracture2D::finalBackboneClustering() {

	N_final_cluster = 0;
	clusterConnection tempcluster;
	tempcluster.LR.first = false;
	tempcluster.LR.second = false;
	for (auto itpci = p_c.begin(); itpci != p_c.end(); itpci++) {
		if ((itpci->cluster_fb == -1)&&(itpci->backboneFinal.size() > 0)) {
			tempcluster.LR.first = false;
			tempcluster.LR.second = false;
			itpci->cluster_fb = ++N_final_cluster;
			if (itpci->left) {
				tempcluster.LR.first = true;
			}
			if (itpci->right) {
				tempcluster.LR.second = true;
			}
			finalCluster.push_back(tempcluster);
			final_clustering(itpci);
		}
	}
}

void fracture2D::final_clustering(list<connectionSPMAT>::iterator &tempPci) {
	

		for (auto jp = tempPci->backboneFinal.begin(); jp != tempPci->backboneFinal.end(); jp++) {
			auto itpcj = p_c.begin();
			advance(itpcj, *jp - 1);
			if (itpcj->cluster_fb == -1) {
				auto itCluster = finalCluster.begin();
				advance(itCluster, tempPci->cluster_fb - 1);
				if ((itpcj->left)&&(!itCluster->LR.first)) {
					itCluster->LR.first = true;
				}
				if ((itpcj->right)&&(!itCluster->LR.second)) {
					itCluster->LR.second = true;
				}
				itpcj->cluster_fb = tempPci->cluster_fb;
				final_clustering(itpcj);
			}
		}

}
			
bool fracture2D::finalConnectedCluster() {
	bool connectivity = false;
	for (auto c = finalCluster.begin(); c != finalCluster.end(); c++) {
		if ((c->LR.first) && (c->LR.second)) {
			connectivity = true;
		}
	}

	for (auto itpci = p_c.begin(); itpci != p_c.end(); itpci++) {
		if ((itpci->cluster_fb > 0) && (itpci->backboneFinal.size() > 0)) {
			auto itCl = finalCluster.begin();
			advance(itCl, itpci->cluster_fb - 1);
			if ((!itCl->LR.first) || (!itCl->LR.second)) {
				while (itpci->backboneFinal.size()){
					itpci->danglingP.push_back(itpci->backboneFinal.front());
					itpci->danglingP_b.push_back(itpci->backboneFinal_b.front());
					itpci->backboneFinal.pop_front();
					itpci->backboneFinal_b.pop_front();
				}
			}
		}
	}

	return connectivity;
}

void fracture2D::printResults() {

	mFstream obj;
	string R_address = txtFolder + "results/" + txtName;
	obj.ResultOutput(R_address);
}

void fracture2D::setResults() {
	results.B = parameter_B;
	results.Bf = parameter_B_final;
	results.cluster1 = NCluster;
	results.cluster2 = int(finalCluster.size());
	results.D = parameter_D;
	results.Df = parameter_Dang;
	results.intersections = int(p_c.size()) - int(initialData.size()) * 2;
	results.K = K1;
	results.Kb = K2;
	results.Lx = Lx;
	results.Ly = Ly;
	results.N = int(initialData.size());
	results.p = occupancy_p;
	results.P = parameter_P;
	results.Pc = Pc;
	results.Q1 = Q_b;
	results.Q2 = Q_b_f;
}

result fracture2D::getResults() {
	return results;
}
