#include "stdafx.h"
#include "CppUnitTest.h"
#include <cmath>
#include"mFstream.h"
#include"global.h"
#include"fracture2D.h"
#include"mException.h"

#include <Eigen/Sparse>




//#include"\Users\Hossein\"


template <typename T>
bool almostEqual(T expected, T actual, T epsilon) {

	return (abs(expected - actual) <= epsilon) ? true : false;
}

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace unitTest
{		
	TEST_CLASS(UnitTest1)
	{
	public:
		
		//TEST_METHOD(TestMeEX)
		//{
		//	double x = 2.1;
		//	double expected = 2.1;
		//	Assert::IsTrue(almostEqual(expected, testMe(x), 1e-6));
		//	// TODO: Your test code here
		//}
		//TEST_METHOD(minEX)
		//{
		//	/*mFstream sample;
		//	char str[MAX_STRING_LENGTH], str1[MAX_STRING_LENGTH];
		//	sample.open("test1.txt", std::fstream::in);

		//	if (!sample.FileSearch("Fractures")); mException("Fractures Not Found");
		//	
		//	if (!sample.FileSearch("Lx")); mException("Lx not found");
		//	if (!sample.ReadWord(str)) mException("Incorrect number for Lx in the data file.");
		//	double Lx = atof(str);

		//	if (!sample.FileSearch("Ly")); mException("Ly not found");
		//	if (!sample.ReadWord(str)) mException("Incorrect number for Ly in the data file.");
		//	double Ly = atof(str);*/



		//	//Assert::IsTrue(almostEqual(expected, minf(x1,x2), 0.0));
		//}

		TEST_METHOD(maxEX)
		{
			double x1 = 10.0;
			double x2 = 2.0;
			double expected = x1;
			Assert::IsTrue(almostEqual(expected, maxf(x1, x2), 0.0));
		}

	};
}


namespace Fracture2D_test
{
	TEST_CLASS(percolation)
	{
	public:

		TEST_METHOD(initialization)
		{
			int N;
			fracture2D percolation2D;

			mFstream sample;
			char str[MAX_STRING_LENGTH], str1[MAX_STRING_LENGTH];
			sample.open("test1.txt", std::fstream::in);
			if (!sample.is_open()) mException("Data file could not be opened.");
			
			if (!sample.FileSearch("NUM")) mException("Number of fractures Not Found");
			if (!sample.ReadWord(str)) mException("Incorrect value for NUM in the data file.");
			N = atoi(str);

			if (!sample.FileSearch("lConst")) mException("lConst not found");
			if (!sample.ReadWord(str)) mException("Incorrect value for lConst in the data file.");
			percolation2D.set_lmin(atof(str));

			if (!sample.FileSearch("Lx")) mException("Lx not found");
			if (!sample.ReadWord(str)) mException("Incorrect value for Lx in the data file.");
			percolation2D.setLx(atof(str));

			if (!sample.FileSearch("Ly")) mException("Ly not found");
			if (!sample.ReadWord(str)) mException("Incorrect value for Ly in the data file.");
			percolation2D.setLy(atof(str));

			if (!sample.FileSearch("Xc")) mException("Xc not found");
			if (!sample.FileSearch("Yc")) mException("Yc not found");
			if (!sample.FileSearch("Orientation")) mException("Orientation not found");

			int i = 0;

			percolation2D.setNcluster(0);
			percolation2D.setTotalPoint(0);

			percolation2D.setConnected(false);
			clusterConnection tempClCon;
			percolation2D.set_a_ex();
			do {
				dataStruct tempInitialData;
				connectionSPMAT tempP_cL, tempP_cR;
				i++;
				percolation2D.setFracNum(i);
				percolation2D.setNcluster(percolation2D.getNcluster() + 1);
				if (!sample.ReadWord(str)) {
					std::cout << "Line : " << N + 4;
					mException("Incorrect value for Xc in the data file.");
				}
				tempInitialData.xCenter = atof(str);

				if (!sample.ReadWord(str)) {
					std::cout << "Line : " << N + 4;
					mException("Incorrect value for Ly in the data file.");
				}
				tempInitialData.yCenter = atof(str);

				if (!sample.ReadWord(str)) {
					std::cout << "Line : " << N + 4;
					mException("Incorrect value for Ly in the data file.");
				}
				tempInitialData.orientation = atof(str);

				cross tempCross;
				tempCross.cluster = percolation2D.getNcluster();
				tempCross.fracNum = i;
				tempCross.connected = false;
				tempClCon.LR.first = false;
				tempClCon.LR.second = false;
				percolation2D.setcluCon(tempClCon);
				percolation2D.setCrossSparseCon(tempCross);
				percolation2D.fractureInitConstantLength(tempInitialData);
				percolation2D.set_initialData(tempInitialData);
				percolation2D.insP_c(tempP_cL, tempP_cR, tempInitialData);
				percolation2D.setp_c(tempP_cL);
				percolation2D.setp_c(tempP_cR);
				percolation2D.findingIntersects();
				percolation2D.calOccupancy_p();
				percolation2D.setConnected(percolation2D.checkLR(tempInitialData));
			} while (i < N);

			bool ans = percolation2D.getConnected();
			Assert::IsTrue(ans);
		}

		TEST_METHOD(test2)
		{
			int N;
			fracture2D percolation2D;

			mFstream sample;
			char str[MAX_STRING_LENGTH], str1[MAX_STRING_LENGTH];
			sample.open("test2.txt", std::fstream::in);
			if (!sample.is_open()) mException("Data file could not be opened.");

			if (!sample.FileSearch("NUM")) mException("Number of fractures Not Found");
			if (!sample.ReadWord(str)) mException("Incorrect value for NUM in the data file.");
			N = atoi(str);

			if (!sample.FileSearch("lConst")) mException("lConst not found");
			if (!sample.ReadWord(str)) mException("Incorrect value for lConst in the data file.");
			percolation2D.set_lmin(atof(str));

			if (!sample.FileSearch("Lx")) mException("Lx not found");
			if (!sample.ReadWord(str)) mException("Incorrect value for Lx in the data file.");
			percolation2D.setLx(atof(str));

			if (!sample.FileSearch("Ly")) mException("Ly not found");
			if (!sample.ReadWord(str)) mException("Incorrect value for Ly in the data file.");
			percolation2D.setLy(atof(str));

			if (!sample.FileSearch("Xc")) mException("Xc not found");
			if (!sample.FileSearch("Yc")) mException("Yc not found");
			if (!sample.FileSearch("Orientation")) mException("Orientation not found");

			int i = 0;

			percolation2D.setNcluster(0);
			percolation2D.setTotalPoint(0);

			percolation2D.setConnected(false);
			clusterConnection tempClCon;
			percolation2D.set_a_ex();
			do {
				dataStruct tempInitialData;
				connectionSPMAT tempP_cL, tempP_cR;
				i++;
				percolation2D.setFracNum(i);
				percolation2D.setNcluster(percolation2D.getNcluster() + 1);
				if (!sample.ReadWord(str)) {
					std::cout << "Line : " << N + 4;
					mException("Incorrect value for Xc in the data file.");
				}
				tempInitialData.xCenter = atof(str);

				if (!sample.ReadWord(str)) {
					std::cout << "Line : " << N + 4;
					mException("Incorrect value for Ly in the data file.");
				}
				tempInitialData.yCenter = atof(str);

				if (!sample.ReadWord(str)) {
					std::cout << "Line : " << N + 4;
					mException("Incorrect value for Orientation in the data file.");
				}
				tempInitialData.orientation = atof(str)*MY_pi/180.0;

				cross tempCross;
				tempCross.cluster = percolation2D.getNcluster();
				tempCross.fracNum = i;
				tempCross.connected = false;
				tempClCon.LR.first = false;
				tempClCon.LR.second = false;
				percolation2D.setcluCon(tempClCon);
				percolation2D.setCrossSparseCon(tempCross);
				percolation2D.fractureInitConstantLength(tempInitialData);
				percolation2D.set_initialData(tempInitialData);
				percolation2D.insP_c(tempP_cL, tempP_cR, tempInitialData);
				percolation2D.setp_c(tempP_cL);
				percolation2D.setp_c(tempP_cR);
				percolation2D.findingIntersects();
				percolation2D.calOccupancy_p();
				percolation2D.setConnected(percolation2D.checkLR(tempInitialData));
			} while (i < N);

			bool ans = percolation2D.getConnected();
			Assert::IsTrue(!ans);
		}


		TEST_METHOD(test3)
		{
			int N;
			fracture2D percolation2D;

			mFstream sample;
			char str[MAX_STRING_LENGTH], str1[MAX_STRING_LENGTH];
			sample.open("test3.txt", std::fstream::in);
			if (!sample.is_open()) mException("Data file could not be opened.");

			if (!sample.FileSearch("NUM")) mException("Number of fractures Not Found");
			if (!sample.ReadWord(str)) mException("Incorrect value for NUM in the data file.");
			N = atoi(str);

			if (!sample.FileSearch("lConst")) mException("lConst not found");
			if (!sample.ReadWord(str)) mException("Incorrect value for lConst in the data file.");
			percolation2D.set_lmin(atof(str));

			if (!sample.FileSearch("Lx")) mException("Lx not found");
			if (!sample.ReadWord(str)) mException("Incorrect value for Lx in the data file.");
			percolation2D.setLx(atof(str));

			if (!sample.FileSearch("Ly")) mException("Ly not found");
			if (!sample.ReadWord(str)) mException("Incorrect value for Ly in the data file.");
			percolation2D.setLy(atof(str));

			if (!sample.FileSearch("Xc")) mException("Xc not found");
			if (!sample.FileSearch("Yc")) mException("Yc not found");
			if (!sample.FileSearch("Orientation")) mException("Orientation not found");

			int i = 0;

			percolation2D.setNcluster(0);
			percolation2D.setTotalPoint(0);

			percolation2D.setConnected(false);
			clusterConnection tempClCon;
			percolation2D.set_a_ex();
			do {
				dataStruct tempInitialData;
				connectionSPMAT tempP_cL, tempP_cR;
				i++;
				percolation2D.setFracNum(i);
				percolation2D.setNcluster(percolation2D.getNcluster() + 1);
				if (!sample.ReadWord(str)) {
					std::cout << "Line : " << N + 4;
					mException("Incorrect value for Xc in the data file.");
				}
				tempInitialData.xCenter = atof(str);

				if (!sample.ReadWord(str)) {
					std::cout << "Line : " << N + 4;
					mException("Incorrect value for Ly in the data file.");
				}
				tempInitialData.yCenter = atof(str);

				if (!sample.ReadWord(str)) {
					std::cout << "Line : " << N + 4;
					mException("Incorrect value for Orientation in the data file.");
				}
				tempInitialData.orientation = atof(str)*MY_pi / 180.0;

				cross tempCross;
				tempCross.cluster = percolation2D.getNcluster();
				tempCross.fracNum = i;
				tempCross.connected = false;
				tempClCon.LR.first = false;
				tempClCon.LR.second = false;
				percolation2D.setcluCon(tempClCon);
				percolation2D.setCrossSparseCon(tempCross);
				percolation2D.fractureInitConstantLength(tempInitialData);
				percolation2D.set_initialData(tempInitialData);
				percolation2D.insP_c(tempP_cL, tempP_cR, tempInitialData);
				percolation2D.setp_c(tempP_cL);
				percolation2D.setp_c(tempP_cR);
				percolation2D.findingIntersects();
				percolation2D.calOccupancy_p();
				percolation2D.setConnected(percolation2D.checkLR(tempInitialData));
			} while (i < N);

			bool ans = percolation2D.getConnected();
			Assert::IsTrue(!ans);
		}

		TEST_METHOD(test4)
		{
			int N;
			fracture2D percolation2D;

			mFstream sample;
			char str[MAX_STRING_LENGTH], str1[MAX_STRING_LENGTH];
			sample.open("test4.txt", std::fstream::in);
			if (!sample.is_open()) mException("Data file could not be opened.");

			if (!sample.FileSearch("NUM")) mException("Number of fractures Not Found");
			if (!sample.ReadWord(str)) mException("Incorrect value for NUM in the data file.");
			N = atoi(str);

			if (!sample.FileSearch("lConst")) mException("lConst not found");
			if (!sample.ReadWord(str)) mException("Incorrect value for lConst in the data file.");
			percolation2D.set_lmin(atof(str));

			if (!sample.FileSearch("Lx")) mException("Lx not found");
			if (!sample.ReadWord(str)) mException("Incorrect value for Lx in the data file.");
			percolation2D.setLx(atof(str));

			if (!sample.FileSearch("Ly")) mException("Ly not found");
			if (!sample.ReadWord(str)) mException("Incorrect value for Ly in the data file.");
			percolation2D.setLy(atof(str));

			if (!sample.FileSearch("Xc")) mException("Xc not found");
			if (!sample.FileSearch("Yc")) mException("Yc not found");
			if (!sample.FileSearch("Orientation")) mException("Orientation not found");

			int i = 0;

			percolation2D.setNcluster(0);
			percolation2D.setTotalPoint(0);

			percolation2D.setConnected(false);
			clusterConnection tempClCon;
			percolation2D.set_a_ex();
			do {
				dataStruct tempInitialData;
				connectionSPMAT tempP_cL, tempP_cR;
				i++;
				percolation2D.setFracNum(i);
				percolation2D.setNcluster(percolation2D.getNcluster() + 1);
				if (!sample.ReadWord(str)) {
					std::cout << "Line : " << N + 4;
					mException("Incorrect value for Xc in the data file.");
				}
				tempInitialData.xCenter = atof(str);

				if (!sample.ReadWord(str)) {
					std::cout << "Line : " << N + 4;
					mException("Incorrect value for Ly in the data file.");
				}
				tempInitialData.yCenter = atof(str);

				if (!sample.ReadWord(str)) {
					std::cout << "Line : " << N + 4;
					mException("Incorrect value for Orientation in the data file.");
				}
				tempInitialData.orientation = atof(str)*MY_pi / 180.0;

				cross tempCross;
				tempCross.cluster = percolation2D.getNcluster();
				tempCross.fracNum = i;
				tempCross.connected = false;
				tempClCon.LR.first = false;
				tempClCon.LR.second = false;
				percolation2D.setcluCon(tempClCon);
				percolation2D.setCrossSparseCon(tempCross);
				percolation2D.fractureInitConstantLength(tempInitialData);
				percolation2D.set_initialData(tempInitialData);
				percolation2D.insP_c(tempP_cL, tempP_cR, tempInitialData);
				percolation2D.setp_c(tempP_cL);
				percolation2D.setp_c(tempP_cR);
				percolation2D.findingIntersects();
				percolation2D.calOccupancy_p();
				percolation2D.setConnected(percolation2D.checkLR(tempInitialData));
			} while (i < N);

			bool ans = percolation2D.getConnected();
			Assert::IsTrue(!ans);
		}

		TEST_METHOD(test5)
		{
			int N;
			fracture2D percolation2D;

			mFstream sample;
			char str[MAX_STRING_LENGTH], str1[MAX_STRING_LENGTH];
			sample.open("test5.txt", std::fstream::in);
			if (!sample.is_open()) mException("Data file could not be opened.");

			if (!sample.FileSearch("NUM")) mException("Number of fractures Not Found");
			if (!sample.ReadWord(str)) mException("Incorrect value for NUM in the data file.");
			N = atoi(str);

			if (!sample.FileSearch("lConst")) mException("lConst not found");
			if (!sample.ReadWord(str)) mException("Incorrect value for lConst in the data file.");
			percolation2D.set_lmin(atof(str));

			if (!sample.FileSearch("Lx")) mException("Lx not found");
			if (!sample.ReadWord(str)) mException("Incorrect value for Lx in the data file.");
			percolation2D.setLx(atof(str));

			if (!sample.FileSearch("Ly")) mException("Ly not found");
			if (!sample.ReadWord(str)) mException("Incorrect value for Ly in the data file.");
			percolation2D.setLy(atof(str));

			if (!sample.FileSearch("Xc")) mException("Xc not found");
			if (!sample.FileSearch("Yc")) mException("Yc not found");
			if (!sample.FileSearch("Orientation")) mException("Orientation not found");

			int i = 0;

			percolation2D.setNcluster(0);
			percolation2D.setTotalPoint(0);

			percolation2D.setConnected(false);
			clusterConnection tempClCon;
			percolation2D.set_a_ex();
			do {
				dataStruct tempInitialData;
				connectionSPMAT tempP_cL, tempP_cR;
				i++;
				percolation2D.setFracNum(i);
				percolation2D.setNcluster(percolation2D.getNcluster() + 1);
				if (!sample.ReadWord(str)) {
					std::cout << "Line : " << N + 4;
					mException("Incorrect value for Xc in the data file.");
				}
				tempInitialData.xCenter = atof(str);

				if (!sample.ReadWord(str)) {
					std::cout << "Line : " << N + 4;
					mException("Incorrect value for Ly in the data file.");
				}
				tempInitialData.yCenter = atof(str);

				if (!sample.ReadWord(str)) {
					std::cout << "Line : " << N + 4;
					mException("Incorrect value for Orientation in the data file.");
				}
				tempInitialData.orientation = atof(str)*MY_pi / 180.0;

				cross tempCross;
				tempCross.cluster = percolation2D.getNcluster();
				tempCross.fracNum = i;
				tempCross.connected = false;
				tempClCon.LR.first = false;
				tempClCon.LR.second = false;
				percolation2D.setcluCon(tempClCon);
				percolation2D.setCrossSparseCon(tempCross);
				percolation2D.fractureInitConstantLength(tempInitialData);
				percolation2D.set_initialData(tempInitialData);
				percolation2D.insP_c(tempP_cL, tempP_cR, tempInitialData);
				percolation2D.setp_c(tempP_cL);
				percolation2D.setp_c(tempP_cR);
				percolation2D.findingIntersects();
				percolation2D.calOccupancy_p();
				percolation2D.setConnected(percolation2D.checkLR(tempInitialData));
			} while (i < N);

			bool ans = percolation2D.getConnected();
			Assert::IsTrue(ans);
		}

		TEST_METHOD(test6)
		{
			int N;
			fracture2D percolation2D;

			mFstream sample;
			char str[MAX_STRING_LENGTH], str1[MAX_STRING_LENGTH];
			sample.open("test6.txt", std::fstream::in);
			if (!sample.is_open()) mException("Data file could not be opened.");

			if (!sample.FileSearch("NUM")) mException("Number of fractures Not Found");
			if (!sample.ReadWord(str)) mException("Incorrect value for NUM in the data file.");
			N = atoi(str);

			if (!sample.FileSearch("lConst")) mException("lConst not found");
			if (!sample.ReadWord(str)) mException("Incorrect value for lConst in the data file.");
			percolation2D.set_lmin(atof(str));

			if (!sample.FileSearch("Lx")) mException("Lx not found");
			if (!sample.ReadWord(str)) mException("Incorrect value for Lx in the data file.");
			percolation2D.setLx(atof(str));

			if (!sample.FileSearch("Ly")) mException("Ly not found");
			if (!sample.ReadWord(str)) mException("Incorrect value for Ly in the data file.");
			percolation2D.setLy(atof(str));

			if (!sample.FileSearch("Xc")) mException("Xc not found");
			if (!sample.FileSearch("Yc")) mException("Yc not found");
			if (!sample.FileSearch("Orientation")) mException("Orientation not found");

			int i = 0;

			percolation2D.setNcluster(0);
			percolation2D.setTotalPoint(0);

			percolation2D.setConnected(false);
			clusterConnection tempClCon;
			percolation2D.set_a_ex();
			do {
				dataStruct tempInitialData;
				connectionSPMAT tempP_cL, tempP_cR;
				i++;
				percolation2D.setFracNum(i);
				percolation2D.setNcluster(percolation2D.getNcluster() + 1);
				if (!sample.ReadWord(str)) {
					std::cout << "Line : " << N + 4;
					mException("Incorrect value for Xc in the data file.");
				}
				tempInitialData.xCenter = atof(str);

				if (!sample.ReadWord(str)) {
					std::cout << "Line : " << N + 4;
					mException("Incorrect value for Ly in the data file.");
				}
				tempInitialData.yCenter = atof(str);

				if (!sample.ReadWord(str)) {
					std::cout << "Line : " << N + 4;
					mException("Incorrect value for Orientation in the data file.");
				}
				tempInitialData.orientation = atof(str)*MY_pi / 180.0;

				cross tempCross;
				tempCross.cluster = percolation2D.getNcluster();
				tempCross.fracNum = i;
				tempCross.connected = false;
				tempClCon.LR.first = false;
				tempClCon.LR.second = false;
				percolation2D.setcluCon(tempClCon);
				percolation2D.setCrossSparseCon(tempCross);
				percolation2D.fractureInitConstantLength(tempInitialData);
				percolation2D.set_initialData(tempInitialData);
				percolation2D.insP_c(tempP_cL, tempP_cR, tempInitialData);
				percolation2D.setp_c(tempP_cL);
				percolation2D.setp_c(tempP_cR);
				percolation2D.findingIntersects();
				percolation2D.calOccupancy_p();
				percolation2D.setConnected(percolation2D.checkLR(tempInitialData));
			} while (i < N);

			bool ans = percolation2D.getConnected();
			Assert::IsTrue(!ans);
		}


		TEST_METHOD(test7)
		{
			int N;
			fracture2D percolation2D;

			mFstream sample;
			char str[MAX_STRING_LENGTH], str1[MAX_STRING_LENGTH];
			sample.open("test7.txt", std::fstream::in);
			if (!sample.is_open()) mException("Data file could not be opened.");

			if (!sample.FileSearch("NUM")) mException("Number of fractures Not Found");
			if (!sample.ReadWord(str)) mException("Incorrect value for NUM in the data file.");
			N = atoi(str);

			if (!sample.FileSearch("lConst")) mException("lConst not found");
			if (!sample.ReadWord(str)) mException("Incorrect value for lConst in the data file.");
			percolation2D.set_lmin(atof(str));

			if (!sample.FileSearch("Lx")) mException("Lx not found");
			if (!sample.ReadWord(str)) mException("Incorrect value for Lx in the data file.");
			percolation2D.setLx(atof(str));

			if (!sample.FileSearch("Ly")) mException("Ly not found");
			if (!sample.ReadWord(str)) mException("Incorrect value for Ly in the data file.");
			percolation2D.setLy(atof(str));

			if (!sample.FileSearch("Xc")) mException("Xc not found");
			if (!sample.FileSearch("Yc")) mException("Yc not found");
			if (!sample.FileSearch("Orientation")) mException("Orientation not found");

			int i = 0;

			percolation2D.setNcluster(0);
			percolation2D.setTotalPoint(0);

			percolation2D.setConnected(false);
			clusterConnection tempClCon;
			percolation2D.set_a_ex();
			do {
				dataStruct tempInitialData;
				connectionSPMAT tempP_cL, tempP_cR;
				i++;
				percolation2D.setFracNum(i);
				percolation2D.setNcluster(percolation2D.getNcluster() + 1);
				if (!sample.ReadWord(str)) {
					std::cout << "Line : " << N + 4;
					mException("Incorrect value for Xc in the data file.");
				}
				tempInitialData.xCenter = atof(str);

				if (!sample.ReadWord(str)) {
					std::cout << "Line : " << N + 4;
					mException("Incorrect value for Ly in the data file.");
				}
				tempInitialData.yCenter = atof(str);

				if (!sample.ReadWord(str)) {
					std::cout << "Line : " << N + 4;
					mException("Incorrect value for Orientation in the data file.");
				}
				tempInitialData.orientation = atof(str)*MY_pi / 180.0;

				cross tempCross;
				tempCross.cluster = percolation2D.getNcluster();
				tempCross.fracNum = i;
				tempCross.connected = false;
				tempClCon.LR.first = false;
				tempClCon.LR.second = false;
				percolation2D.setcluCon(tempClCon);
				percolation2D.setCrossSparseCon(tempCross);
				percolation2D.fractureInitConstantLength(tempInitialData);
				percolation2D.set_initialData(tempInitialData);
				percolation2D.insP_c(tempP_cL, tempP_cR, tempInitialData);
				percolation2D.setp_c(tempP_cL);
				percolation2D.setp_c(tempP_cR);
				percolation2D.findingIntersects();
				percolation2D.calOccupancy_p();
				percolation2D.setConnected(percolation2D.checkLR(tempInitialData));
			} while (i < N);
			percolation2D.updateFlags();
			percolation2D.updateCross();
			percolation2D.RunPercolation();
			percolation2D.percolationResult();
			bool ans = percolation2D.getConnected();
			Assert::IsTrue(ans);


		}


		TEST_METHOD(test8)
		{
			int N;
			fracture2D percolation2D;

			mFstream sample;
			char str[MAX_STRING_LENGTH], str1[MAX_STRING_LENGTH];
			sample.open("test8.txt", std::fstream::in);
			if (!sample.is_open()) mException("Data file could not be opened.");

			if (!sample.FileSearch("NUM")) mException("Number of fractures Not Found");
			if (!sample.ReadWord(str)) mException("Incorrect value for NUM in the data file.");
			N = atoi(str);

			if (!sample.FileSearch("lConst")) mException("lConst not found");
			if (!sample.ReadWord(str)) mException("Incorrect value for lConst in the data file.");
			percolation2D.set_lmin(atof(str));

			if (!sample.FileSearch("Lx")) mException("Lx not found");
			if (!sample.ReadWord(str)) mException("Incorrect value for Lx in the data file.");
			percolation2D.setLx(atof(str));

			if (!sample.FileSearch("Ly")) mException("Ly not found");
			if (!sample.ReadWord(str)) mException("Incorrect value for Ly in the data file.");
			percolation2D.setLy(atof(str));

			if (!sample.FileSearch("Xc")) mException("Xc not found");
			if (!sample.FileSearch("Yc")) mException("Yc not found");
			if (!sample.FileSearch("Orientation")) mException("Orientation not found");

			int i = 0;

			percolation2D.setNcluster(0);
			percolation2D.setTotalPoint(0);

			percolation2D.setConnected(false);
			clusterConnection tempClCon;
			percolation2D.set_a_ex();
			do {
				dataStruct tempInitialData;
				connectionSPMAT tempP_cL, tempP_cR;
				i++;
				percolation2D.setFracNum(i);
				percolation2D.setNcluster(percolation2D.getNcluster() + 1);
				if (!sample.ReadWord(str)) {
					std::cout << "Line : " << N + 4;
					mException("Incorrect value for Xc in the data file.");
				}
				tempInitialData.xCenter = atof(str);

				if (!sample.ReadWord(str)) {
					std::cout << "Line : " << N + 4;
					mException("Incorrect value for Ly in the data file.");
				}
				tempInitialData.yCenter = atof(str);

				if (!sample.ReadWord(str)) {
					std::cout << "Line : " << N + 4;
					mException("Incorrect value for Orientation in the data file.");
				}
				tempInitialData.orientation = atof(str)*MY_pi / 180.0;

				cross tempCross;
				tempCross.cluster = percolation2D.getNcluster();
				tempCross.fracNum = i;
				tempCross.connected = false;
				tempClCon.LR.first = false;
				tempClCon.LR.second = false;
				percolation2D.setcluCon(tempClCon);
				percolation2D.setCrossSparseCon(tempCross);
				percolation2D.fractureInitConstantLength(tempInitialData);
				percolation2D.set_initialData(tempInitialData);
				percolation2D.insP_c(tempP_cL, tempP_cR, tempInitialData);
				percolation2D.setp_c(tempP_cL);
				percolation2D.setp_c(tempP_cR);
				percolation2D.findingIntersects();
				percolation2D.calOccupancy_p();
				percolation2D.setConnected(percolation2D.checkLR(tempInitialData));
			} while (i < N);
			percolation2D.updateFlags();
			percolation2D.updateCross();
			percolation2D.RunPercolation();
			percolation2D.percolationResult();
			percolation2D.networkResult();
			bool ans = percolation2D.getConnected();
			Assert::IsTrue(ans);


		}


		TEST_METHOD(test9)
		{
			int N;
			fracture2D percolation2D;

			mFstream sample;
			char str[MAX_STRING_LENGTH], str1[MAX_STRING_LENGTH];
			sample.open("test9.txt", std::fstream::in);
			if (!sample.is_open()) mException("Data file could not be opened.");

			if (!sample.FileSearch("NUM")) mException("Number of fractures Not Found");
			if (!sample.ReadWord(str)) mException("Incorrect value for NUM in the data file.");
			N = atoi(str);

			if (!sample.FileSearch("lConst")) mException("lConst not found");
			if (!sample.ReadWord(str)) mException("Incorrect value for lConst in the data file.");
			percolation2D.set_lmin(atof(str));

			if (!sample.FileSearch("Lx")) mException("Lx not found");
			if (!sample.ReadWord(str)) mException("Incorrect value for Lx in the data file.");
			percolation2D.setLx(atof(str));

			if (!sample.FileSearch("Ly")) mException("Ly not found");
			if (!sample.ReadWord(str)) mException("Incorrect value for Ly in the data file.");
			percolation2D.setLy(atof(str));

			if (!sample.FileSearch("Xc")) mException("Xc not found");
			if (!sample.FileSearch("Yc")) mException("Yc not found");
			if (!sample.FileSearch("Orientation")) mException("Orientation not found");

			int i = 0;

			percolation2D.setNcluster(0);
			percolation2D.setTotalPoint(0);

			percolation2D.setConnected(false);
			clusterConnection tempClCon;
			percolation2D.set_a_ex();
			do {
				dataStruct tempInitialData;
				connectionSPMAT tempP_cL, tempP_cR;
				i++;
				percolation2D.setFracNum(i);
				percolation2D.setNcluster(percolation2D.getNcluster() + 1);
				if (!sample.ReadWord(str)) {
					std::cout << "Line : " << N + 4;
					mException("Incorrect value for Xc in the data file.");
				}
				tempInitialData.xCenter = atof(str);

				if (!sample.ReadWord(str)) {
					std::cout << "Line : " << N + 4;
					mException("Incorrect value for Ly in the data file.");
				}
				tempInitialData.yCenter = atof(str);

				if (!sample.ReadWord(str)) {
					std::cout << "Line : " << N + 4;
					mException("Incorrect value for Orientation in the data file.");
				}
				tempInitialData.orientation = atof(str)*MY_pi / 180.0;

				cross tempCross;
				tempCross.cluster = percolation2D.getNcluster();
				tempCross.fracNum = i;
				tempCross.connected = false;
				tempClCon.LR.first = false;
				tempClCon.LR.second = false;
				percolation2D.setcluCon(tempClCon);
				percolation2D.setCrossSparseCon(tempCross);
				percolation2D.fractureInitConstantLength(tempInitialData);
				percolation2D.set_initialData(tempInitialData);
				percolation2D.insP_c(tempP_cL, tempP_cR, tempInitialData);
				percolation2D.setp_c(tempP_cL);
				percolation2D.setp_c(tempP_cR);
				percolation2D.findingIntersects();
				percolation2D.calOccupancy_p();
				percolation2D.setConnected(percolation2D.checkLR(tempInitialData));
			} while (i < N);
			percolation2D.updateFlags();
			percolation2D.updateCross();
			percolation2D.RunPercolation();
			percolation2D.percolationResult();
			percolation2D.networkResult();
			bool ans = percolation2D.getConnected();
			Assert::IsTrue(ans);


		}

		TEST_METHOD(test10)
		{
			int N;
			fracture2D percolation2D;

			mFstream sample;
			char str[MAX_STRING_LENGTH], str1[MAX_STRING_LENGTH];
			sample.open("test10.txt", std::fstream::in);
			if (!sample.is_open()) mException("Data file could not be opened.");

			if (!sample.FileSearch("NUM")) mException("Number of fractures Not Found");
			if (!sample.ReadWord(str)) mException("Incorrect value for NUM in the data file.");
			N = atoi(str);

			if (!sample.FileSearch("lConst")) mException("lConst not found");
			if (!sample.ReadWord(str)) mException("Incorrect value for lConst in the data file.");
			percolation2D.set_lmin(atof(str));

			if (!sample.FileSearch("Lx")) mException("Lx not found");
			if (!sample.ReadWord(str)) mException("Incorrect value for Lx in the data file.");
			percolation2D.setLx(atof(str));

			if (!sample.FileSearch("Ly")) mException("Ly not found");
			if (!sample.ReadWord(str)) mException("Incorrect value for Ly in the data file.");
			percolation2D.setLy(atof(str));

			if (!sample.FileSearch("Xc")) mException("Xc not found");
			if (!sample.FileSearch("Yc")) mException("Yc not found");
			if (!sample.FileSearch("Orientation")) mException("Orientation not found");

			int i = 0;

			percolation2D.setNcluster(0);
			percolation2D.setTotalPoint(0);

			percolation2D.setConnected(false);
			clusterConnection tempClCon;
			percolation2D.set_a_ex();
			do {
				dataStruct tempInitialData;
				connectionSPMAT tempP_cL, tempP_cR;
				i++;
				percolation2D.setFracNum(i);
				percolation2D.setNcluster(percolation2D.getNcluster() + 1);
				if (!sample.ReadWord(str)) {
					std::cout << "Line : " << N + 4;
					mException("Incorrect value for Xc in the data file.");
				}
				tempInitialData.xCenter = atof(str);

				if (!sample.ReadWord(str)) {
					std::cout << "Line : " << N + 4;
					mException("Incorrect value for Ly in the data file.");
				}
				tempInitialData.yCenter = atof(str);

				if (!sample.ReadWord(str)) {
					std::cout << "Line : " << N + 4;
					mException("Incorrect value for Orientation in the data file.");
				}
				tempInitialData.orientation = atof(str)*MY_pi / 180.0;

				cross tempCross;
				tempCross.cluster = percolation2D.getNcluster();
				tempCross.fracNum = i;
				tempCross.connected = false;
				tempClCon.LR.first = false;
				tempClCon.LR.second = false;
				percolation2D.setcluCon(tempClCon);
				percolation2D.setCrossSparseCon(tempCross);
				percolation2D.fractureInitConstantLength(tempInitialData);
				percolation2D.set_initialData(tempInitialData);
				percolation2D.insP_c(tempP_cL, tempP_cR, tempInitialData);
				percolation2D.setp_c(tempP_cL);
				percolation2D.setp_c(tempP_cR);
				percolation2D.findingIntersects();
				percolation2D.calOccupancy_p();
				percolation2D.setConnected(percolation2D.checkLR(tempInitialData));
			} while (i < N);
			percolation2D.updateFlags();
			percolation2D.updateCross();
			percolation2D.RunPercolation();
			percolation2D.percolationResult();
			percolation2D.networkResult();
			bool ans = percolation2D.getConnected();
			Assert::IsTrue(ans);


		}


		TEST_METHOD(test11)
		{
			int N;
			fracture2D percolation2D;

			mFstream sample;
			char str[MAX_STRING_LENGTH], str1[MAX_STRING_LENGTH];
			sample.open("test11.txt", std::fstream::in);
			if (!sample.is_open()) mException("Data file could not be opened.");

			if (!sample.FileSearch("NUM")) mException("Number of fractures Not Found");
			if (!sample.ReadWord(str)) mException("Incorrect value for NUM in the data file.");
			N = atoi(str);

			if (!sample.FileSearch("lConst")) mException("lConst not found");
			if (!sample.ReadWord(str)) mException("Incorrect value for lConst in the data file.");
			percolation2D.set_lmin(atof(str));

			if (!sample.FileSearch("Lx")) mException("Lx not found");
			if (!sample.ReadWord(str)) mException("Incorrect value for Lx in the data file.");
			percolation2D.setLx(atof(str));

			if (!sample.FileSearch("Ly")) mException("Ly not found");
			if (!sample.ReadWord(str)) mException("Incorrect value for Ly in the data file.");
			percolation2D.setLy(atof(str));

			if (!sample.FileSearch("Xc")) mException("Xc not found");
			if (!sample.FileSearch("Yc")) mException("Yc not found");
			if (!sample.FileSearch("Orientation")) mException("Orientation not found");

			int i = 0;

			percolation2D.setNcluster(0);
			percolation2D.setTotalPoint(0);

			percolation2D.setConnected(false);
			clusterConnection tempClCon;
			percolation2D.set_a_ex();
			do {
				dataStruct tempInitialData;
				connectionSPMAT tempP_cL, tempP_cR;
				i++;
				percolation2D.setFracNum(i);
				percolation2D.setNcluster(percolation2D.getNcluster() + 1);
				if (!sample.ReadWord(str)) {
					std::cout << "Line : " << N + 4;
					mException("Incorrect value for Xc in the data file.");
				}
				tempInitialData.xCenter = atof(str);

				if (!sample.ReadWord(str)) {
					std::cout << "Line : " << N + 4;
					mException("Incorrect value for Ly in the data file.");
				}
				tempInitialData.yCenter = atof(str);

				if (!sample.ReadWord(str)) {
					std::cout << "Line : " << N + 4;
					mException("Incorrect value for Orientation in the data file.");
				}
				tempInitialData.orientation = atof(str)*MY_pi / 180.0;

				cross tempCross;
				tempCross.cluster = percolation2D.getNcluster();
				tempCross.fracNum = i;
				tempCross.connected = false;
				tempClCon.LR.first = false;
				tempClCon.LR.second = false;
				percolation2D.setcluCon(tempClCon);
				percolation2D.setCrossSparseCon(tempCross);
				percolation2D.fractureInitConstantLength(tempInitialData);
				percolation2D.set_initialData(tempInitialData);
				percolation2D.insP_c(tempP_cL, tempP_cR, tempInitialData);
				percolation2D.setp_c(tempP_cL);
				percolation2D.setp_c(tempP_cR);
				percolation2D.findingIntersects();
				percolation2D.calOccupancy_p();
				percolation2D.setConnected(percolation2D.checkLR(tempInitialData));
			} while (i < N);
			percolation2D.updateFlags();
			percolation2D.updateCross();
			percolation2D.RunPercolation();
			percolation2D.percolationResult();
			percolation2D.networkResult();
			bool ans = percolation2D.getConnected();
			Assert::IsTrue(ans);


		}


	};
}

