#pragma once
#define tol_q 0.02
#include<chrono>
#include<random>
#include<cmath>
#include<list>
#include<algorithm>
#include<iostream>
//#include <string>
#include<chrono>
#include"global.h"
#include"mFstream.h"
#include"mException.h"
#include <Eigen/Sparse>
#include<Eigen/IterativeLinearSolvers>
#include <Eigen/Dense>

using Eigen::SparseMatrix; 
using namespace Eigen;

using namespace std;

class fracture2D
{
public:
	
	fracture2D(int n);
	fracture2D(double p, double lmin, double lmax, double a, double tethaMin, double tethaMax, double Lx, double Ly, string folder, string fileName, double Pr, double Pl);
	fracture2D(double p, double lf, double tethaMin, double tethaMax, double Lx, double Ly, string folder, string fileName, double Pr, double Pl);
	fracture2D();



	~fracture2D();

	double lengthPowerLawDist();
	double fractureOrientation();
	double fractureCenterLocation(double L);
	double fractureDelY(double orientation, double length);
	double fractureDelX(double orientation, double length);
	void get_initialData(list<dataStruct> &inputStruc);
	double calLength(double x1, double x2, double y1, double y2);
	void fractureInit(dataStruct & temp);
	void fractureInitConstantLength(dataStruct & temp);
	void limitBoundary(dataStruct &temp);
	void calAveExArea(list<dataStruct>::iterator lastElement);
	double calExArea(double l);
	void calOccupancy_p();
	void sysEffLength();
	void updateCSPMAT(connectionSPMAT insData);
	void findingIntersects();
	bool interstionValidation(pair<double, double> fL1, pair<double, double> fR1, pair<double, double> fL2, pair<double, double> fR2, pair<double, double> P, int i, int j);
	//void mergingOverlapLines(list<dataStruct>::iterator beginning, list<dataStruct>::iterator ending);
	double uniformDist();
	void insP_c(connectionSPMAT &outL, connectionSPMAT &outR, dataStruct &in);
	list<connectionSPMAT> get_p_c(list<connectionSPMAT> in_p_c) const;
	list<connectionSPMAT> set_p_c(list<connectionSPMAT> out_p_c) ;
	bool checkLR(const dataStruct &tempInitialData);
	void mergePoints();
	void RunPercolation();
	bool p_c_equality(connectionSPMAT p_c_1, connectionSPMAT p_c_2);
	connectionSPMAT merge_p_c(connectionSPMAT &temp1, connectionSPMAT &temp2);
	void deleteingIsolatedFractures();
	void deleteingIsolatedFractures_d();
	//void mergeSubCon(subCon &subCon1, subCon &subCon2);
	void findNearestPoint();
	void percolationProbability();
	void backboneProbability();
	void backboneFinalProbability();
	void deletingDeadEndFractures();
	void deletingDeadEndFractures_d();
	void findingDeadends(list<connectionSPMAT>::iterator &tempPci);
	void findingDeadends_d(list<connectionSPMAT>::iterator &tempPci);
	void updateFlags();
	void setLx(double Lx);
	void setLy(double Ly);
	void set_lmin(double lmin);
	void set_a(double a);
	void set_initialData(dataStruct temp);
	void set_p_c(connectionSPMAT tempP_c);
	void setTotalPoint(int n);
	void setNcluster(int n);
	void setConnected(bool val);
	double calAbsPerm();
	double calAbsPerm_final();
	int getTotalPoint();
	int getNcluster();
	void randomFractureGen(dataStruct &temp);
	bool getConnected();
	void setcluCon(clusterConnection &temp);
	void setCrossSparseCon(cross &temp);
	void setp_c(connectionSPMAT &temp);
	void set_a_ex();
	void setFracNum(int N);
	void mergeLine(int line1, int line2);
	bool checkMergeLine(int line1, int line2);
	void percolationResult();
	void updateCross();
	void networkResult();
	void backboneResult();
	void backbone_finalResult();
	void backbone_indexing();
	void backbone_indexing_final();
	void pressureAssignment_backbone(VectorXd x);
	void deleting_dangling();
	void finalBackboneClustering();
	void final_clustering(list<connectionSPMAT>::iterator &tempPci);
	bool finalConnectedCluster();
	void printResults();
	void setResults();
	result getResults();

private:
	//fracture2D();
	double Q_b;
	double Q_b_f;
	int N;
	double p;
	double lmin;
	double lmax;
	double a;
	double tethaMin;
	double tethaMax;
	double Lx;
	double Ly;
	double a_ex;
	double occupancy_p;
	double effLx;
	double effLy;
	double aveProjlX;
	int totalPointNum;
	double aveProjlY;
	int NCluster;
	bool connected;
	list<dataStruct> initialData;
	dataStruct tempInitialData;
	list<connectionSPMAT> p_c;
	//list<connectionSPMAT> p_c_merged;
	list<cross> crossSparse;
	list<clusterConnection> cluCon;
	double parameter_P;
	double parameter_B;
	double parameter_D;
	double parameter_B_final;
	double parameter_Dang;
	double S_OP_l;
	double S_PerP_l;
	double S_PerB_l;
	double S_PerD_l;
	double S_PerB_lf;
	double S_PerD_lf;
	string txtFolder;
	string txtName;
	int nnz;
	int nnz_final;

	double PR, PL;
	int numFrBac;
	int numFrBac_final;
	list<clusterConnection> finalCluster;
	int N_final_cluster;

	double K1, K2;
	result results;
	int intersections;
	double Pc;

};

