#pragma once

#include<list>

using namespace std;

#define MY_pi 3.141592653589793
#define equality_Tol 1e-32
#define DarcyToM2 0.9866*1e-12
#define mToft 3.28084

struct result {
	double Lx;
	double Ly;
	double N;
	double p;
	double P;
	double B;
	double D;
	double Bf;
	double Df;
	double K;
	double Kb;
	int iteration;
	int cluster1;
	int cluster2;
	int intersections;
	double Pc;
	double Q1;
	double Q2;
};



struct dataStruct {
	double xCenter;
	double yCenter;
	double orientation;
	double b;
	double length;
	double xR;
	double yR;
	double xL;
	double yL;
	double slope; // y = slope*x + c;
	double c;
	pair<int, int> p_n;
};

struct subCon
{
	list<int> pointJ;
	double slope;
	double c;
	double b;
};

struct connectionSPMAT {
	int pointNum;
	list<subCon> conecction;
	double xi;
	double yi;
	//bool merged;
	int cluster;
	int cluster_fb;
	bool left;
	bool right;
	bool flag;
	bool deadend;
	bool dangling;
	//bool merge;
	list<int> nearestPoint;
	list<int> backbone;
	list<int> backboneFinal;
	list<int> deadEndP;
	list<int> danglingP;

	list<double> nearestPoint_b;
	list<double> backbone_b;
	list<double> backboneFinal_b;
	list<double> deadEndP_b;
	list<double> danglingP_b;

	int i_b;
	int i_b_f;

	double P_b;
	double P_b_f;

};

struct cross
{
	int fracNum;
	list<int> col;
	int cluster;
	bool connected;
};

struct clusterConnection
{
	pair<bool, bool> LR;
};


static double minf(double x1, double x2) {
	return x1 < x2 ? x1 : x2;
}

static double almostEqual(double &x1, double &x2, double tol) {

	if (abs((x1 - x2)/x1) < tol) {
		return true;
	}
	else {
		return false;
	}
}

static double maxf(double x1, double x2) {
	return x1 > x2 ? x1 : x2;
}

