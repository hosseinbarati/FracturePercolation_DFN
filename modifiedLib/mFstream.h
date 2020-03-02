#pragma once
#include<fstream>
#include<iomanip>
#include<string>
#include<list>
#include"global.h"
#include <sstream>
#include "mException.h"
#include <iostream>


#define MAX_STRING_LENGTH 256

using namespace std;

class mFstream : public fstream 
{
public:
	mFstream();
	~mFstream();
	mFstream(result results, string address);
	void writeInitaialState(list<dataStruct> initialStat);
	void write_p_c(connectionSPMAT it, int N);
	void writeCluster(list<cross> clusterList, list<dataStruct> initialData);
	void onlineWriteCluster(list<cross> clusterList, int N, list<dataStruct> initialData);
	void exceptionsOutput(string msg);
	int FileSearch(char* );
	int ReadWord(char* );
	void percolatedOutput(list<cross> &crossSparse, list<dataStruct> &initialData, string address);
	void networkOutput(list<dataStruct> &initialData, string address);
	void backboneOutput(list<connectionSPMAT> &p_c, string address);
	void backbonefinalOutput(list<connectionSPMAT> &p_c, string address);
	void ResultOutput(string address);
	void printResults(result results, string address, int fileNum, int N);
	
private:
	static int lineintWrite;
	ofstream ofile;
	ofstream clusterText;

	ofstream Exceptions;
	ofstream summary;

};

