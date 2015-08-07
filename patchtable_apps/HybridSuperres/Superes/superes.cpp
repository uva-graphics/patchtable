// SupeRes.cpp : Defines the entry point for the console application.TrainImgColsize
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <stdio.h>      /* printf */
#include <assert.h>     /* assert */
#include<stdlib.h>
#include "time.h"
#include "patchtable.h"
using namespace cv;
#include <opencv2/imgproc/imgproc.hpp>
int PCADIMS = 10;
using namespace std;
string TRAINFILENAME;
vector<string> TESTFILENAME;
string OUTPUTFILENAME;
Mat TrainImg, TrainCR, TrainCB, TrainLum, *TestImg, *TestLum, *TestCR, *TestCB, *ResultImg;
Mat ** HighRSPatch;
int vectorSize = 256;
int TrainImgHRSize = 72, TrainImgLRSize = 8, DownSampleFactor = 9;
Mat OriginalTrainFeatureData;
PCA *pca;
Mat TrainFeatureData;
int KeyRows = 0;
int TrainImgColsize = 0;
int TrainImgRowsize = 0;
//#define assert2(_Expression, _Msg) (void)( (!!(_Expression)) || (_wassert(_CRT_WIDE(#_Msg), _CRT_WIDE(__FILE__), __LINE__), 0) )
#define assert2(expr, msg) ASSERT2(expr, msg)
string TEST = "images\\";
//#define DrawResult 1
//#define BigPicture 1
#define Patchtable 1
#define SmallPicture 1
void WriteFile(int rows, int cols, float **data, string NAME){
	float x_min = 500.0, x_max = -500.0;
	ofstream out(NAME);
	for (int i = 0; i < rows; i++){
		for (int j = 0; j < cols - 1; j++){
			out << data[i][j] << " ";
			if (data[i][j] > x_max)
			{
				x_max = data[i][j];
			}
			if (data[i][j] < x_min)
			{
				x_min = data[i][j];
			}
		}
		out << data[i][cols - 1] << endl;
		if (data[i][cols - 1] > x_max)
		{
			x_max = data[i][cols - 1];
		}
		if (data[i][cols - 1] < x_min)
		{
			x_min = data[i][cols - 1];
		}
	}
	out << x_min;
	out << x_max;
	out.close();
}
void LoadImg(string testname, string trainame, string outputname)
{
	int positionPotAfter = testname.find_first_of(".");
	string testPotAfter = testname.substr(0, positionPotAfter);
	string testForMat = testname.substr(positionPotAfter + 1, testname.size());
	if (testForMat == "txt")
	{
		ifstream in(testname);
		char line[1024];
		string nameTemp;
		stringstream ss2;
		while (in.getline(line, 1024)){
			ss2 << line;
			ss2 >> nameTemp;
			TESTFILENAME.push_back(nameTemp);
		}
		in.close();
	}
	else
	{
		TESTFILENAME.push_back(testname);
	}
	OUTPUTFILENAME = outputname;
	int testNo = TESTFILENAME.size();
	TestImg = new  Mat[testNo];
	TestCB = new  Mat[testNo];
	TestCR = new  Mat[testNo];
	TestLum = new  Mat[testNo];
	vector< Mat> LumYcYb;
	for (int i = 0; i < testNo; i++)
	{
		TestImg[i] = imread(TESTFILENAME[i], 1);
		TestImg[i].convertTo(TestImg[i], CV_32FC3);
		TestImg[i] = TestImg[i] / 255.0;
		cvtColor(TestImg[i], TestImg[i], COLOR_BGR2YCrCb);
		split(TestImg[i], LumYcYb);
		TestLum[i] = LumYcYb[0].clone();
		TestCR[i] = LumYcYb[1].clone();
		TestCB[i] = LumYcYb[2].clone();
		LumYcYb.clear();
	}
	TrainImg = imread(TRAINFILENAME, 1);
	TrainImg.convertTo(TrainImg, CV_32FC3);
	TrainImg = TrainImg / 255.0;
	cvtColor(TrainImg, TrainImg, COLOR_BGR2YCrCb);
	split(TrainImg, LumYcYb);
	TrainLum = LumYcYb[0].clone();
	TrainCR = LumYcYb[1].clone();
	TrainCB = LumYcYb[2].clone();
	LumYcYb.clear();
	delete[]TestImg;
}
void Is_boundary(int &row, int &col, int rowpathsize, int colpathsize, int filter)
{
	if (filter == 1 || filter == 3)
	{
		col = col < 0 ? 0 : col;
		col = (col > colpathsize - 1) ? (colpathsize - 1) : col;
	}
	else if (filter == 2 || filter == 4)
	{
		row = (row > rowpathsize - 1) ? (rowpathsize - 1) : row;
		row = row < 0 ? 0 : row;
	}
}
void Filter(Mat lowpatch1, Mat & random_vector, int mode)
{
	int row1 = 0, col1 = 0, row2 = 0, col2 = 0, n = 0; float a = 0.0, b = 0.0;
	if (mode == 1)
	{
		for (int r = 0; r < lowpatch1.rows; r++)
		{
			for (int c = 0; c < lowpatch1.cols; c++)
			{
				row1 = r, col1 = c - 1, row2 = r, col2 = c + 1;
				Is_boundary(row1, col1, lowpatch1.rows, lowpatch1.cols, 1);
				a = lowpatch1.at<float>(row1, col1);
				Is_boundary(row2, col2, lowpatch1.rows, lowpatch1.cols, 1);
				b = lowpatch1.at<float>(row2, col2);
				random_vector.at<float>(r, c) = b - a;
			}
		}
	}
	else if (mode == 2)
	{
		for (int r = 0; r < lowpatch1.rows; r++)
		{
			for (int c = 0; c < lowpatch1.cols; c++)
			{
				row1 = r - 1, col1 = c, row2 = r + 1, col2 = c;
				Is_boundary(row1, col1, lowpatch1.rows, lowpatch1.cols, 2);
				a = lowpatch1.at<float>(row1, col1);
				Is_boundary(row2, col2, lowpatch1.rows, lowpatch1.cols, 2);
				b = lowpatch1.at<float>(row2, col2);
				random_vector.at<float>(r, c) = b - a;
			}
		}
	}
	else if (mode == 3)
	{
		for (int r = 0; r < lowpatch1.rows; r++)
		{
			for (int c = 0; c < lowpatch1.cols; c++)
			{
				row1 = r, col1 = c - 2, row2 = r, col2 = c + 2;
				Is_boundary(row1, col1, lowpatch1.rows, lowpatch1.cols, mode);
				a = lowpatch1.at<float>(row1, col1);
				Is_boundary(row2, col2, lowpatch1.rows, lowpatch1.cols, mode);
				b = lowpatch1.at<float>(row2, col2);
				random_vector.at<float>(r, c) = b + a - 2 * lowpatch1.at<float>(r, c);
			}
		}
	}
	else if (mode == 4)
	{
		for (int r = 0; r < lowpatch1.rows; r++)
		{
			for (int c = 0; c < lowpatch1.cols; c++)
			{
				row1 = r - 2, col1 = c, row2 = r + 2, col2 = c;
				Is_boundary(row1, col1, lowpatch1.rows, lowpatch1.cols, mode);
				a = lowpatch1.at<float>(row1, col1);
				Is_boundary(row2, col2, lowpatch1.rows, lowpatch1.cols, mode);
				b = lowpatch1.at<float>(row2, col2);
				random_vector.at<float>(r, c) = b + a - 2 * lowpatch1.at<float>(r, c);
			}
		}
	}
}
void Concatenat(Mat lowpatch1, Mat & random_vector, int keyrow)
{
	int row1 = 0, col1 = 0, row2 = 0, col2 = 0, n = 0; float a = 0.0, b = 0.0;
	for (int r = 0; r < TrainImgLRSize; r++)
	{
		for (int c = 0; c < TrainImgLRSize; c++)
		{
			row1 = r, col1 = c - 1, row2 = r, col2 = c + 1;
			Is_boundary(row1, col1, TrainImgLRSize, TrainImgLRSize, 1);
			a = lowpatch1.at<float>(row1, col1);
			Is_boundary(row2, col2, TrainImgLRSize, TrainImgLRSize, 1);
			b = lowpatch1.at<float>(row2, col2);
			random_vector.at<float>(keyrow, n) = b - a; n++;
		}
	}
	for (int r = 0; r < TrainImgLRSize; r++)
	{
		for (int c = 0; c < TrainImgLRSize; c++)
		{
			row1 = r - 1, col1 = c, row2 = r + 1, col2 = c;
			Is_boundary(row1, col1, TrainImgLRSize, TrainImgLRSize, 2);
			a = lowpatch1.at<float>(row1, col1);
			Is_boundary(row2, col2, TrainImgLRSize, TrainImgLRSize, 2);
			b = lowpatch1.at<float>(row2, col2);
			random_vector.at<float>(keyrow, n) = b - a; n++;
		}
	}
	for (int r = 0; r < TrainImgLRSize; r++)
	{
		for (int c = 0; c < TrainImgLRSize; c++)
		{
			row1 = r, col1 = c - 2, row2 = r, col2 = c + 2;
			Is_boundary(row1, col1, TrainImgLRSize, TrainImgLRSize, 3);
			a = lowpatch1.at<float>(row1, col1);
			Is_boundary(row2, col2, TrainImgLRSize, TrainImgLRSize, 3);
			b = lowpatch1.at<float>(row2, col2);
			random_vector.at<float>(keyrow, n) = b + a - 2 * lowpatch1.at<float>(r, c); n++;
		}
	}
	for (int r = 0; r < TrainImgLRSize; r++)
	{
		for (int c = 0; c < TrainImgLRSize; c++)
		{
			row1 = r - 2, col1 = c, row2 = r + 2, col2 = c;
			Is_boundary(row1, col1, TrainImgLRSize, TrainImgLRSize, 4);
			a = lowpatch1.at<float>(row1, col1);
			Is_boundary(row2, col2, TrainImgLRSize, TrainImgLRSize, 4);
			b = lowpatch1.at<float>(row2, col2);
			random_vector.at<float>(keyrow, n) = b + a - 2 * lowpatch1.at<float>(r, c); n++;
		}
	}
}
//just back up
//float* read_points(const char* filename, int rows, int cols)
//{
//	float* data;
//	float *p;
//	FILE* fin;
//	int i, j;
//
//	fin = fopen(filename, "r");
//	if (!fin) {
//		printf("Cannot open input file.\n");
//		exit(1);
//	}
//
//	data = (float*)malloc(rows*cols*sizeof(float));
//	if (!data) {
//		printf("Cannot allocate memory.\n");
//		exit(1);
//	}
//	p = data;
//
//	for (i = 0; i < rows; ++i) {
//		for (j = 0; j < cols; ++j) {
//			fscanf(fin, "%f ", p);
//			p++;
//		}
//	}
//	fclose(fin);
//	return data;
//}
template<typename T>
void write_results(const char* filename, T **data, int rows, int cols)
{
	ofstream out(filename);
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			if (typeid(T).name() == typeid(1.0f).name())
			{
				out << data[i][j] << endl;
			}
			else if (typeid(T).name() == typeid(1).name())
			{
				out << data[i][j] << endl;
			}
			else
			{
				assert2(typeid(T).name() == typeid(1.0f).name() || typeid(T).name() == typeid(1).name(), "expected int or float");
			}
		}
	}
	out.close();
}
template<typename T>
void ReaDate(int rows, int cols, T** dataset, const char* filename)
{
	ifstream in(filename);
	char line[200];
	int lineno = 0; int i = 0, j = 0;
	while (in.getline(line, 200)){
		stringstream ss2;
		ss2 << line;
		ss2 >> dataset[i][j];
		j++;
		if (j == cols)
		{
			j = 0; i++;
		}
	}
	in.close();
}
#ifdef DrawResult
float TrainMaxValue = -1000, TrainMinValue = 1000;
#endif
void GeneraTrainSet()
{
	//int TrainImgRowsize = 0, TrainImgColsize = 0;
	TrainImgColsize = TrainLum.cols;
	TrainImgRowsize = TrainLum.rows;
	KeyRows = (TrainImgColsize - TrainImgHRSize + 1)*(TrainImgRowsize - TrainImgHRSize + 1);//highres patch number
#ifdef BigPicture
	Mat randomKey(10000, vectorSize, CV_32FC1);
	//imwrite("TrainImg.png", TrainLum*255.0);
	/*Adapt to big data*/
	RNG rng;
	int i = 0;
	Mat temp = Mat::zeros(TrainImgRowsize - TrainImgHRSize + 1, TrainImgColsize - TrainImgHRSize + 1, CV_32SC1);
	do{
		int col = rng.uniform(0, TrainImgColsize - TrainImgHRSize + 1);
		int row = rng.uniform(0, TrainImgRowsize - TrainImgHRSize + 1);
		if (temp.at<int>(row, col) == 0)
		{
			temp.at<int>(row, col) = 1;
			Mat highpatch(TrainImgHRSize, TrainImgHRSize, CV_32FC1);
			Mat lowpatch1 =  Mat::zeros(TrainImgLRSize, TrainImgLRSize, CV_32FC1);
			TrainLum( Rect(col, row, TrainImgHRSize, TrainImgHRSize)).copyTo(highpatch);
			resize(highpatch, lowpatch1, Size(highpatch.cols / DownSampleFactor, highpatch.rows / DownSampleFactor), 0, 0, CV_INTER_CUBIC);
			Concatenat(lowpatch1, randomKey,i);
			i++;
		}
	} while (i < 10000);
	pca = new  PCA(randomKey, Mat(), CV_PCA_DATA_AS_ROW, PCADIMS);
	TrainFeatureData = Mat::zeros(KeyRows, PCADIMS, CV_32FC1);
#else
	OriginalTrainFeatureData = Mat::zeros(KeyRows, vectorSize, CV_32FC1);
#endif
#pragma omp parallel for
	for (int r = 0; r < TrainImgRowsize - TrainImgHRSize + 1; r++)
	{
		for (int c = 0; c < TrainImgColsize - TrainImgHRSize + 1; c++)
		{
			Mat input1(1, vectorSize, CV_32FC1);
			Mat highpatch(TrainImgHRSize, TrainImgHRSize, CV_32FC1);
			Mat lowpatch1 = Mat::zeros(TrainImgLRSize, TrainImgLRSize, CV_32FC1);
			TrainLum(Rect(c, r, TrainImgHRSize, TrainImgHRSize)).copyTo(highpatch);
			resize(highpatch, lowpatch1, Size(highpatch.cols / DownSampleFactor, highpatch.rows / DownSampleFactor), 0, 0, INTER_AREA);
			Concatenat(lowpatch1, input1, 0);
#ifdef BigPicture
			Mat input2(1, PCADIMS, CV_32FC1);
			pca->project(input1, input2);
			for (int n = 0; n < input2.cols; n++)
			{
				TrainFeatureData.at<float>(r*(TrainImgColsize - TrainImgHRSize + 1) + c, n) = input2.at<float>(0, n);
#ifdef DrawResult
				TrainMaxValue = (input1.at<float>(0, n) - TrainMaxValue)>0 ? input1.at<float>(0, n) : TrainMaxValue;
				TrainMinValue = (TrainMinValue - input1.at<float>(0, n))>0 ? input1.at<float>(0, n) : TrainMinValue;
#endif
			}
#else
			for (int n = 0; n < input1.cols; n++)
			{
				OriginalTrainFeatureData.at<float>(r*(TrainImgColsize - TrainImgHRSize + 1) + c, n) = input1.at<float>(0, n);
#ifdef DrawResult
				TrainMaxValue = (input1.at<float>(0, n) - TrainMaxValue)>0 ? input1.at<float>(0, n) : TrainMaxValue;
				TrainMinValue = (TrainMinValue - input1.at<float>(0, n))>0 ? input1.at<float>(0, n) : TrainMinValue;
#endif
			}

#endif
		}
	}
#ifndef BigPicture
	pca = new  PCA(OriginalTrainFeatureData, Mat(), CV_PCA_DATA_AS_ROW, PCADIMS);
	TrainFeatureData = Mat::zeros(KeyRows, PCADIMS, CV_32FC1);
	pca->project(OriginalTrainFeatureData, TrainFeatureData);
#endif

}
void InWhichPic(int index, int lowRows, int &r, int &c)
{
	int n = 0;
	for (int x = 0; x < lowRows; x++)
	{
		if (index >= lowRows*x && index < lowRows*(x + 1))
		{
			r = x;
			c = index - x*lowRows;
			break;
		}

	}

}
float Average(Mat& input)
{
	float sum = 0.0;
	for (int r = 0; r < input.rows; r++)
	{
		for (int c = 0; c < input.cols; c++)
		{
			sum = sum + input.at<float>(r, c);
		}
	}
	sum = sum / (input.rows*input.cols);
	return sum;
}
void SuperRe()
{
	int nn = 9;
	int** result;
	int testNo = TESTFILENAME.size();
	ResultImg = new  Mat[testNo];
	int index = 0, resultIndex = 0;
	clock_t start, finish;
	double  duration;
	start = clock();
#ifndef Patchtable
	//////////////////////////////FOR OPENCV FLANN
	flann::Index flann_index(
		TrainFeatureData,
		flann::KDTreeIndexParams(4),
		cvflann::FLANN_DIST_L2
		);
	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	cout << duration << ".seconds"<<endl;
	vector<int> indices(nn);
	vector<float> dists(nn);
#endif

#ifdef DrawResult
	float MaxValue = -10000.0, MinValue = 100000.0;
#endif
	for (int itestNo = 0; itestNo < testNo; itestNo++)
	{
		int lowRows = TestLum[itestNo].rows;
		int lowCols = TestLum[itestNo].cols;
		int tcount = (lowRows - TrainImgLRSize + 1)*(lowCols - TrainImgLRSize + 1);
		result = new int*[tcount];
		for (int n = 0; n < tcount; n++)
		{
			result[n] = new int[nn];
		}
		Mat filter1mat = Mat::zeros(lowRows, lowCols, CV_32FC1);
		Mat filter2mat = Mat::zeros(lowRows, lowCols, CV_32FC1);
		Mat filter3mat = Mat::zeros(lowRows, lowCols, CV_32FC1);
		Mat filter4mat = Mat::zeros(lowRows, lowCols, CV_32FC1);
		Filter(TestLum[itestNo], filter1mat, 1);
		Filter(TestLum[itestNo], filter2mat, 2);
		Filter(TestLum[itestNo], filter3mat, 3);
		Filter(TestLum[itestNo], filter4mat, 4);
		ResultImg[itestNo] = Mat::zeros(lowRows*DownSampleFactor, lowCols*DownSampleFactor, CV_32FC3);
		Mat finalresult = Mat::zeros(lowRows*DownSampleFactor, lowCols*DownSampleFactor, CV_32FC3);
		Mat SplitTest = Mat::zeros((lowRows - TrainImgLRSize + 1)*TrainImgLRSize, (lowCols - TrainImgLRSize + 1)*TrainImgLRSize, CV_32FC1);
		Mat testemp(TrainImgLRSize, TrainImgLRSize, CV_32FC1);
		int coutcount = 0;
#ifdef Patchtable
		Mat TestFeatureData = Mat::zeros(tcount, PCADIMS, CV_32FC1);
#endif
#ifdef DrawResult
		for (int r = 0; r < lowRows - TrainImgLRSize + 1; r++)
		{
			for (int c = 0; c < lowCols - TrainImgLRSize + 1; c++)
			{
				Mat input1(1, vectorSize, CV_32FC1);
				Mat output1(1, 20, CV_32FC1);
				Mat testMatch = Mat::ones(TrainImgHRSize + 1, (8 + 1 + 8 + 1 + 8 + 1 + 8 + 72 + 1), CV_32FC1);
				TestLum[itestNo](Rect(c, r, TrainImgLRSize, TrainImgLRSize)).copyTo(testemp);
				testemp.copyTo(SplitTest(Rect(c*TrainImgLRSize, r*TrainImgLRSize, TrainImgLRSize, TrainImgLRSize)));
				filter1mat(Rect(c, r, TrainImgLRSize, TrainImgLRSize)).copyTo(testemp);
				int n = 0;
				for (int i = 0; i < TrainImgLRSize; i++)
				{
					for (int j = 0; j < TrainImgLRSize; j++)
					{
						input1.at<float>(0, n) = testemp.at<float>(i, j);
						if (input1.at<float>(0, n) > MaxValue)
						{
							MaxValue = input1.at<float>(0, n);
						}
						if (input1.at<float>(0, n) < MinValue)
						{
							MinValue = input1.at<float>(0, n);
						}
						n++;
					}
				}
				filter2mat(Rect(c, r, TrainImgLRSize, TrainImgLRSize)).copyTo(testemp);
				for (int i = 0; i < TrainImgLRSize; i++)
				{
					for (int j = 0; j < TrainImgLRSize; j++)
					{
						input1.at<float>(0, n) = testemp.at<float>(i, j);
						if (input1.at<float>(0, n) > MaxValue)
						{
							MaxValue = input1.at<float>(0, n);
						}
						if (input1.at<float>(0, n) < MinValue)
						{
							MinValue = input1.at<float>(0, n);
						}
						n++;
					}
				}
				filter3mat(Rect(c, r, TrainImgLRSize, TrainImgLRSize)).copyTo(testemp);
				for (int i = 0; i < TrainImgLRSize; i++)
				{
					for (int j = 0; j < TrainImgLRSize; j++)
					{
						input1.at<float>(0, n) = testemp.at<float>(i, j);
						if (input1.at<float>(0, n) > MaxValue)
						{
							MaxValue = input1.at<float>(0, n);
						}
						if (input1.at<float>(0, n) < MinValue)
						{
							MinValue = input1.at<float>(0, n);
						}
						n++;
					}
				}
				filter4mat(Rect(c, r, TrainImgLRSize, TrainImgLRSize)).copyTo(testemp);
				for (int i = 0; i < TrainImgLRSize; i++)
				{
					for (int j = 0; j < TrainImgLRSize; j++)
					{
						input1.at<float>(0, n) = testemp.at<float>(i, j);
						if (input1.at<float>(0, n) > MaxValue)
						{
							MaxValue = input1.at<float>(0, n);
						}
						if (input1.at<float>(0, n) < MinValue)
						{
							MinValue = input1.at<float>(0, n);
						}
						n++;
					}
				}
			}
		}
#endif
		duration = 0.0;
		for (int r = 0; r < lowRows - TrainImgLRSize + 1; r++)
		{
			for (int c = 0; c < lowCols - TrainImgLRSize + 1; c++)
			{
				//r = 6, c = 11;
				Mat input1(1, vectorSize, CV_32FC1);
				Mat output1(1, PCADIMS, CV_32FC1);
				Mat testMatch = Mat::ones(TrainImgHRSize + 1, (8 + 1 + 8 + 1 + 8 + 1 + 8 + 72 + 1), CV_32FC1);
				TestLum[itestNo](Rect(c, r, TrainImgLRSize, TrainImgLRSize)).copyTo(testemp);
				testemp.copyTo(SplitTest(Rect(c*TrainImgLRSize, r*TrainImgLRSize, TrainImgLRSize, TrainImgLRSize)));
				filter1mat(Rect(c, r, TrainImgLRSize, TrainImgLRSize)).copyTo(testemp);
				int n = 0;
				for (int i = 0; i < TrainImgLRSize; i++)
				{
					for (int j = 0; j < TrainImgLRSize; j++)
					{
						//input1.at<float>(0, n) = (testemp.at<float>(i, j) - MinValue) / (MaxValue - MinValue);
						input1.at<float>(0, n) = testemp.at<float>(i, j);
						n++;
					}
				}
				filter2mat(Rect(c, r, TrainImgLRSize, TrainImgLRSize)).copyTo(testemp);
				for (int i = 0; i < TrainImgLRSize; i++)
				{
					for (int j = 0; j < TrainImgLRSize; j++)
					{
						//input1.at<float>(0, n) = (testemp.at<float>(i, j) - MinValue) / (MaxValue - MinValue);
						input1.at<float>(0, n) = testemp.at<float>(i, j);
						n++;
					}
				}
				filter3mat(Rect(c, r, TrainImgLRSize, TrainImgLRSize)).copyTo(testemp);
				for (int i = 0; i < TrainImgLRSize; i++)
				{
					for (int j = 0; j < TrainImgLRSize; j++)
					{
						//input1.at<float>(0, n) = (testemp.at<float>(i, j) - MinValue) / (MaxValue - MinValue);
						input1.at<float>(0, n) = testemp.at<float>(i, j);
						n++;
					}
				}
				filter4mat(Rect(c, r, TrainImgLRSize, TrainImgLRSize)).copyTo(testemp);
				for (int i = 0; i < TrainImgLRSize; i++)
				{
					for (int j = 0; j < TrainImgLRSize; j++)
					{
						//input1.at<float>(0, n) = (testemp.at<float>(i, j) - MinValue) / (MaxValue - MinValue);
						input1.at<float>(0, n) = testemp.at<float>(i, j);
						n++;
					}
				}
				pca->project(input1, output1);
#ifdef Patchtable
				for (int q = 0; q < output1.cols; q++)
				{
					TestFeatureData.at<float>(coutcount, q) = output1.at<float>(0, q);
				}
#else
				start = clock();
				flann_index.knnSearch(output1, indices, dists, nn, flann::SearchParams(256));
				finish = clock();
				duration = (double)(finish - start) / CLOCKS_PER_SEC + duration;
				vector<int>::iterator indexn = indices.begin(); int x2 = 0; int matchr = 0, matchc = 0;
				vector<float>::iterator indexndis = dists.begin();
				float disMin = 10000.0;
				for (; indexn != indices.end(); indexn++)
				{
					int n = *indexn;
#ifdef DrawResult
					InWhichPic(n, TrainLum.cols - TrainImgHRSize + 1, matchr, matchc);
					int testr = 0, testc = 0;// float distance = 0.0;
					for (int j = 0; j < vectorSize; j++) {
						testMatch.at<float>(testr, testc) = (input1.at<float>(0, j) - MinValue) / (MaxValue - MinValue);
						testMatch.at<float>(testr, testc + 9) = (OriginalTrainFeatureData.at<float>(matchr*(TrainImgColsize - TrainImgHRSize + 1) + matchc, j) - TrainMinValue) / (TrainMaxValue - TrainMinValue);
						//cout << testMatch.at<float>(testr, testc) << " " << testMatch.at<float>(testr, testc + 9) << " " << testMatch.at<float>(testr, testc) - testMatch.at<float>(testr, testc + 9) << endl;
						//distance = abs(testMatch.at<float>(testr, testc) - testMatch.at<float>(testr, testc + 9)) + distance; testc++;
						if (testc == 8)
						{
							testc = 0; testr++;
						}
					}
					//cout << "the distance of the original is" << 0.88 << ";   the distance of the our result is" << distance;
					for (int tesr = 0; tesr < TrainImgLRSize; tesr++)
					{
						for (int tesc = 0; tesc < TrainImgLRSize; tesc++)
						{
							testMatch.at<float>(tesr, 18 + tesc) = TestLum[itestNo].at<float>(r + tesr, c + tesc);
						}
					}
					Mat testemphr(TrainImgHRSize, TrainImgHRSize, CV_32FC1);
					for (int r = 0; r < TrainImgHRSize; r++)
					{
						for (int c = 0; c < TrainImgHRSize; c++)
						{
							testMatch.at<float>(0 + r, 36 + c) = TrainLum.at<float>(matchr + r, matchc + c);
							testemphr.at<float>(r, c) = TrainLum.at<float>(matchr + r, matchc + c);
						}
					}
					resize(testemphr, testemp, Size(TrainImgLRSize, TrainImgLRSize), 0, 0, INTER_CUBIC);
					for (int r = 0; r < TrainImgLRSize; r++)
					{
						for (int c = 0; c < TrainImgLRSize; c++)
						{
							testMatch.at<float>(0 + r, 27 + c) = testemp.at<float>(r, c);
						}
					}
					char file_name[15];
					sprintf(file_name, "TestFirstMatch%d.png", n);
					cv::imwrite(file_name, testMatch*255.0);
#endif
					if (*indexndis < disMin)
					{
						disMin = *indexndis;
						result[coutcount][x2] = result[coutcount][0]; x2++;
						result[coutcount][0] = n;
					}
					else
					{
						result[coutcount][x2] = n; x2++;
					}
					indexndis++;
				}
				indices.clear();
				dists.clear();
#endif
				coutcount++;
			}
		}
#ifndef Patchtable
		cout<<"search time:"<<duration<<".seconds";
#else
		//dealwith TrainFeatureData

		int trainh = TrainImgRowsize - TrainImgHRSize + 1;
		int trainw = TrainImgColsize - TrainImgHRSize + 1;
		Array<float> TrainVectorArray(trainh, trainw, PCADIMS);

		int count = 0;
		for (int r = 0; r < trainh; r++)
		{
			for (int c = 0; c < trainw; c++)
			{
				for (int k = 0; k < PCADIMS; k++)
				{
					TrainVectorArray(r, c, k) = TrainFeatureData.at<float>(count,k);
				}
				count++;
			}
		}
		int testh = lowRows - TrainImgLRSize + 1;
		int testw = lowCols - TrainImgLRSize + 1;
		Array<float> TestVectorArray(testh, testw, PCADIMS);
		count = 0;
		for (int r = 0; r < testh; r++)
		{
			for (int c = 0; c < testw; c++)
			{
				for (int k = 0; k < PCADIMS; k++)
				{
					TestVectorArray(r, c, k) = TestFeatureData.at<float>(count, k);
				}
				count++;
			}
		}
		PatchTableParams *p = new PatchTableParams();
		p->set_speed(0);                    // Needs to be above the other lines
		//p->coherence_spatial = 4.0;
		//p->calc_exact_dist = true;
		p->ndims = PCADIMS;
		p->nchroma = 0;
		p->is_descriptor = true;

		//p->limit = 4e6;    // Number of cells in table
		//p->do_rs = true;   // Whether to run random search
		//p->run_dt = true;  // Whether to fill in empty cells with distance transform
		//p->do_prop = true; // Whether to run propagation

		//allow contains the TrainFeatureData which is height*width rows and PCADIMS cols

		PatchTable<>* table = new PatchTable<>(p, TrainVectorArray);
		Array<double> ann;
		
		double latest_time = 1e100;
		for (int iter = 0; iter < 2; iter++){
			latest_time = table->lookup(TestVectorArray, ann);
		}

#endif
		//ReaDate(KeyRows, vectorSize, testset, "testset.txt");
		//ReaDate(tcount, nn, result, "result.txt");
		Mat weightResult = Mat::zeros(lowRows*DownSampleFactor, lowCols*DownSampleFactor, CV_32SC1);
		Mat Result = Mat::zeros(lowRows*DownSampleFactor, lowCols*DownSampleFactor, CV_32FC1);
		int r = 0, c = 0, tempR = 0, tempC = 0; float sumweight = 0.0; int index[9];
		Mat resultFindHighPatch = Mat::zeros(TrainImgHRSize * 50, TrainImgHRSize*(nn), CV_32FC1);
		Mat resultFindLowPatch = Mat::zeros(TrainImgLRSize * 50, TrainImgLRSize*(nn + 1), CV_32FC1);

		int ttcount = 0;
		for (int r = 0; r < (lowRows - TrainImgLRSize + 1); r++)
		{
			for (int c = 0; c < (lowCols - TrainImgLRSize + 1); c++)
			{
				//r = 6, c = 11;
				Mat tempHighPatchSum = Mat::zeros(TrainImgHRSize, TrainImgHRSize, CV_32FC1);
				Mat tempHighPatch = Mat::zeros(TrainImgHRSize, TrainImgHRSize, CV_32FC1);
				Mat tempLowPatch = Mat::zeros(TrainImgLRSize, TrainImgLRSize, CV_32FC1);
				TestLum[0](Range(r, r + TrainImgLRSize), Range(c, c + TrainImgLRSize)).copyTo(tempLowPatch);
				for (int i = 0; i < 1; i++)
				{
#ifndef	Patchtable
					index[i] = result[ttcount][i];
					InWhichPic(index[i], (TrainImg.cols - TrainImgHRSize + 1), tempR, tempC);
#else
					tempC = ann(r, c, NNF_X);
					tempR = ann(r, c, NNF_Y);
#ifdef DrawResult
					Mat input1(1, vectorSize, CV_32FC1);
					filter1mat(Rect(c, r, TrainImgLRSize, TrainImgLRSize)).copyTo(testemp);
					int n = 0;
					for (int i = 0; i < TrainImgLRSize; i++)
					{
						for (int j = 0; j < TrainImgLRSize; j++)
						{
							//input1.at<float>(0, n) = (testemp.at<float>(i, j) - MinValue) / (MaxValue - MinValue);
							input1.at<float>(0, n) = testemp.at<float>(i, j);
							n++;
						}
					}
					filter2mat(Rect(c, r, TrainImgLRSize, TrainImgLRSize)).copyTo(testemp);
					for (int i = 0; i < TrainImgLRSize; i++)
					{
						for (int j = 0; j < TrainImgLRSize; j++)
						{
							//input1.at<float>(0, n) = (testemp.at<float>(i, j) - MinValue) / (MaxValue - MinValue);
							input1.at<float>(0, n) = testemp.at<float>(i, j);
							n++;
						}
					}
					filter3mat(Rect(c, r, TrainImgLRSize, TrainImgLRSize)).copyTo(testemp);
					for (int i = 0; i < TrainImgLRSize; i++)
					{
						for (int j = 0; j < TrainImgLRSize; j++)
						{
							//input1.at<float>(0, n) = (testemp.at<float>(i, j) - MinValue) / (MaxValue - MinValue);
							input1.at<float>(0, n) = testemp.at<float>(i, j);
							n++;
						}
					}
					filter4mat(Rect(c, r, TrainImgLRSize, TrainImgLRSize)).copyTo(testemp);
					for (int i = 0; i < TrainImgLRSize; i++)
					{
						for (int j = 0; j < TrainImgLRSize; j++)
						{
							//input1.at<float>(0, n) = (testemp.at<float>(i, j) - MinValue) / (MaxValue - MinValue);
							input1.at<float>(0, n) = testemp.at<float>(i, j);
							n++;
						}
					}
					Mat testMatch = Mat::ones(TrainImgHRSize + 1, (8 + 1 + 8 + 1 + 8 + 1 + 8 + 72 + 1), CV_32FC1);
					TestLum[itestNo](Rect(c, r, TrainImgLRSize, TrainImgLRSize)).copyTo(testemp);
					testemp.copyTo(SplitTest(Rect(c*TrainImgLRSize, r*TrainImgLRSize, TrainImgLRSize, TrainImgLRSize)));
					int testr = 0, testc = 0; //float distance = 0;
					for (int j = 0; j < vectorSize; j++) {
						testMatch.at<float>(testr, testc) = (input1.at<float>(0, j) - MinValue) / (MaxValue - MinValue);
						testMatch.at<float>(testr, testc + 9) = (OriginalTrainFeatureData.at<float>(tempR*(TrainImgColsize - TrainImgHRSize + 1) + tempC, j) - TrainMinValue) / (TrainMaxValue - TrainMinValue);
						//cout << testMatch.at<float>(testr, testc) << " " << testMatch.at<float>(testr, testc + 9) << " " << testMatch.at<float>(testr, testc) - testMatch.at<float>(testr, testc + 9)<<endl;
						//distance = abs(testMatch.at<float>(testr, testc) - testMatch.at<float>(testr, testc + 9))+distance;
						testc++;
						if (testc == 8)
						{
							testc = 0; testr++;
						}
					}
					//cout << "the distance of the original is" << 13.4538 << ";   the distance of the our result is" << distance;
					for (int tesr = 0; tesr < TrainImgLRSize; tesr++)
					{
						for (int tesc = 0; tesc < TrainImgLRSize; tesc++)
						{
							testMatch.at<float>(tesr, 18 + tesc) = TestLum[itestNo].at<float>(r + tesr, c + tesc);
						}
					}
					Mat testemphr(TrainImgHRSize, TrainImgHRSize, CV_32FC1);
					for (int r = 0; r < TrainImgHRSize; r++)
					{
						for (int c = 0; c < TrainImgHRSize; c++)
						{
							testMatch.at<float>(0 + r, 36 + c) = TrainLum.at<float>(tempR + r, tempC + c);
							testemphr.at<float>(r, c) = TrainLum.at<float>(tempR + r, tempC + c);
						}
					}
					resize(testemphr, testemp, Size(TrainImgLRSize, TrainImgLRSize), 0, 0, INTER_CUBIC);
					for (int r = 0; r < TrainImgLRSize; r++)
					{
						for (int c = 0; c < TrainImgLRSize; c++)
						{
							testMatch.at<float>(0 + r, 27 + c) = testemp.at<float>(r, c);
						}
					}
					char file_name[15];
					sprintf(file_name, "TestFirstMatch%d.png", tempR*(TrainImgColsize - TrainImgHRSize + 1) + tempC);
					cv::imwrite(file_name, testMatch*255.0);
#endif
#endif
					TrainLum(Range(tempR, tempR + TrainImgHRSize), Range(tempC, tempC + TrainImgHRSize)).copyTo(tempHighPatch);
					float average_highpatch = Average(tempHighPatch);
					float average_lowpatch = Average(tempLowPatch);
					//resize(tempHighPatch, tempLowPatch, Size(tempHighPatch.cols / DownSampleFactor, tempHighPatch.rows / DownSampleFactor), 0, 0, CV_INTER_CUBIC);
					//tempHighPatchSum = exp(param) * tempHighPatch + tempHighPatchSum;
					tempHighPatchSum = tempHighPatch - average_highpatch + average_lowpatch;
				}
				//tempHighPatchSum = tempHighPatchSum / sumweight;
				/*char file_name[15];
				sprintf(file_name, "TestFinal%d.png", ttcount);
				imwrite(file_name, tempHighPatchSum*255.0);*/
				ttcount++;
				sumweight = 0, 0;
				int startIndexR = r * DownSampleFactor; int startIndexC = c * DownSampleFactor;
				for (int i = 0; i < TrainImgHRSize; i++)
				{
					for (int j = 0; j < TrainImgHRSize; j++)
					{
						Result.at<float>(i + startIndexR, j + startIndexC) = tempHighPatchSum.at<float>(i, j) + Result.at<float>(i + startIndexR, j + startIndexC);
						weightResult.at<int>(i + startIndexR, j + startIndexC) = 1 + weightResult.at<int>(i + startIndexR, j + startIndexC);
					}
				}
			}
		}
		for (int i = 0; i < Result.rows; i++)
		{
			for (int j = 0; j < Result.cols; j++)
			{
				Result.at<float>(i, j) = Result.at<float>(i, j) / weightResult.at<int>(i, j);
			}
		}
		imwrite(TEST + "addOriginalResult.png", Result*255.0);
		TestLum[itestNo] = Result;
		Mat CR = Result.clone();
		Mat CB = Result.clone();
		Mat CRt = Mat::zeros(CR.rows, CR.cols, CV_32FC1);
		Mat CBt = Mat::zeros(CR.rows, CR.cols, CV_32FC1);
		resize(TestCR[itestNo], CR, Size(TestCB[itestNo].cols * DownSampleFactor, TestCB[itestNo].rows * DownSampleFactor), 0, 0, INTER_CUBIC);
		bilateralFilter(CR, CRt, -1.0, 0.1, DownSampleFactor + 2);
		resize(TestCB[itestNo], CB, Size(TestCB[itestNo].cols * DownSampleFactor, TestCB[itestNo].rows * DownSampleFactor), 0, 0, INTER_CUBIC);
		bilateralFilter(CB, CBt, -1.0, 0.1, DownSampleFactor + 2);
		vector< Mat> LumYcYb;
		LumYcYb.push_back(Result.clone());
		LumYcYb.push_back(CRt.clone());
		LumYcYb.push_back(CBt.clone());
		merge(LumYcYb, ResultImg[itestNo]);
		cvtColor(ResultImg[itestNo], ResultImg[itestNo], COLOR_YCrCb2BGR);
		ResultImg[itestNo] = ResultImg[itestNo] * 255.0;
		imwrite(TEST + "addfinalResult.jpg", ResultImg[itestNo]);
		for (int m = 0; m < tcount; m++)
		{
			delete[] result[m];
		}
		delete[] result;
	}
	delete[]TestLum;
	delete[]ResultImg;
	delete[]TestImg;
	delete[]TestCB;
	delete[]TestCR;
//	Sleep(5000);
}
void usage()
{
	printf("Superes highResImg.png view_01_01.png output.png\n");
	printf("Or if you have a list names of png in a txt file, then you can enter：\n");
	printf("Superes highResImg.png view.txt output.png\n");
}
void Test1()
{
	int TrainImgRowsize = 0, TrainImgColsize = 0;
	TrainImgColsize = TrainLum.cols;
	TrainImgRowsize = TrainLum.rows;
	Mat tempSplitPatches = Mat::zeros(((TrainImgRowsize - TrainImgHRSize) / DownSampleFactor + 1)*TrainImgHRSize, ((TrainImgColsize - TrainImgHRSize) / DownSampleFactor + 1)*TrainImgHRSize, CV_32FC1);
	Mat tempMatchpatches = Mat::zeros(((TrainImgRowsize - TrainImgHRSize) / DownSampleFactor + 1)*((TrainImgColsize - TrainImgHRSize) / DownSampleFactor + 1)*TrainImgLRSize, 2 * TrainImgLRSize, CV_32FC1);
	int i = 0, j = 0, indexj = 0, indexi = 0, testi = 0, testj = 0;
	for (int r = 0; r < TrainImgRowsize - TrainImgHRSize; r = r + DownSampleFactor)
	{
		for (int c = 0; c < TrainImgColsize - TrainImgHRSize; c = c + DownSampleFactor)
		{
			Mat highpatch(TrainImgHRSize, TrainImgHRSize, CV_32FC1);
			Mat lowpatch1(TrainImgLRSize, TrainImgLRSize, CV_32FC1);
			Mat FirstPatch(TrainImgLRSize, TrainImgLRSize, CV_32FC1);
			TrainLum(Rect(c, r, TrainImgHRSize, TrainImgHRSize)).copyTo(highpatch);
			if (r == 21 * DownSampleFactor&&c == 17 * DownSampleFactor)
			{
				int reeoe = 0;
				imwrite("TargetHRImg.png", highpatch*255.0);
			}
			highpatch.copyTo(tempSplitPatches(Rect(j*TrainImgHRSize, i*TrainImgHRSize, TrainImgHRSize, TrainImgHRSize)));
			resize(highpatch, lowpatch1, Size(highpatch.cols / DownSampleFactor, highpatch.rows / DownSampleFactor), 0, 0, INTER_AREA);
			if (r == 21 * DownSampleFactor&&c == 17 * DownSampleFactor)
			{
				int reeoe = 0;
				imwrite("TargetLRImg.png", lowpatch1*255.0);
			}
			for (int p = 0; p < TrainImgLRSize; p++)
			{
				for (int q = 0; q < TrainImgLRSize; q++)
				{
					float eroor = lowpatch1.at<float>(p, q);
					float yy = eroor;
				}
			}
			Mat input1(1, vectorSize, CV_32FC1);
			Mat input2(1, vectorSize, CV_32FC1);
			Concatenat(lowpatch1, input1, 0);
			for (int p = 0; p < vectorSize; p++)
			{
				float eroor = input1.at<float>(0, p);
				float yy = eroor;
			}
			TestLum[0](Rect(testj, testi, TrainImgLRSize, TrainImgLRSize)).copyTo(FirstPatch);
			if (r == 21 * DownSampleFactor&&c == 17 * DownSampleFactor)
			{
				int reeoe = 0;
				imwrite("TestLRImg.png", FirstPatch*255.0);
			}
			for (int p = 0; p < TrainImgLRSize; p++)
			{
				for (int q = 0; q < TrainImgLRSize; q++)
				{
					float eroor = FirstPatch.at<float>(p, q);
					float yy = eroor;
				}
			}
			Concatenat(FirstPatch, input2, 0);
			for (int p = 0; p < vectorSize; p++)
			{
				float eroor = input1.at<float>(0, p);
				float yy = eroor;
			}
			float delta = 0.0;
			for (int p = 0; p < vectorSize; p++)
			{
				delta += (input1.at<float>(0, p) - input2.at<float>(0, p))*(input1.at<float>(0, p) - input2.at<float>(0, p));
				float yy = delta;
			}
			float yy = delta;
			TestLum[0](Rect(testj, testi, TrainImgLRSize, TrainImgLRSize)).copyTo(tempMatchpatches(Rect(indexj*TrainImgLRSize, indexi*TrainImgLRSize, TrainImgLRSize, TrainImgLRSize)));
			indexj++;
			lowpatch1.copyTo(tempMatchpatches(Rect(indexj*TrainImgLRSize, indexi*TrainImgLRSize, TrainImgLRSize, TrainImgLRSize)));
			j++;
			indexi++;
			indexj = 0;
			testj++;
		}
		i++; j = 0; testj = 0; testi++;
	}
	imwrite("FunALLHighPatch.png", tempSplitPatches*255.0);
	imwrite("FunFindLowHighPatchMatch.png", tempMatchpatches*255.0);
	int lowRows = TestLum[0].rows;
	int lowCols = TestLum[0].cols;
	Mat weightResult = Mat::zeros(lowRows*DownSampleFactor, lowCols*DownSampleFactor, CV_32SC1);
	Mat Result = Mat::zeros(lowRows*DownSampleFactor, lowCols*DownSampleFactor, CV_32FC1);
	i = 0, j = 0;
	for (int r = 0; r < Result.rows - TrainImgHRSize; r = r + DownSampleFactor)
	{
		for (int c = 0; c < Result.cols - TrainImgHRSize; c = c + DownSampleFactor)
		{
			for (int indexi = 0; indexi < TrainImgHRSize; indexi++)
			{
				for (int indexj = 0; indexj < TrainImgHRSize; indexj++)
				{
					Result.at<float>(r + indexi, c + indexj) += tempSplitPatches.at<float>(i + indexi, j + indexj);
					weightResult.at<int>(r + indexi, c + indexj) += 1;

				}
			}
			j = j + TrainImgHRSize;
		}
		i = i + TrainImgHRSize; j = 0;
	}
	for (int i = 0; i < Result.rows; i++)
	{
		for (int j = 0; j < Result.cols; j++)
		{
			Result.at<float>(i, j) = Result.at<float>(i, j) / weightResult.at<int>(i, j);
		}
	}
	ofstream inlum("my_lum.txt");
	for (int i = 0; i < Result.rows; i++)
	{
		for (int j = 0; j < Result.cols; j++)
		{
			inlum << Result.at<float>(i, j) << " ";
		}
		inlum << endl;
	}
	inlum.close();
	imwrite(TEST + "addOriginalResult.png", Result*255.0);
}
void Test2()
{
	Mat test = Mat::zeros(TrainImgHRSize, TrainImgHRSize, CV_32FC1);
	TrainLum(Rect(396, 102, TrainImgHRSize, TrainImgHRSize)).copyTo(test);
	imwrite("64656_1.png", test*255.0);
	TrainLum(Rect(396, 103, TrainImgHRSize, TrainImgHRSize)).copyTo(test);
	imwrite("64656_2.png", test*255.0);
	TrainLum(Rect(396, 101, TrainImgHRSize, TrainImgHRSize)).copyTo(test);
	imwrite("64656_3.png", test*255.0);
	TrainLum(Rect(397, 102, TrainImgHRSize, TrainImgHRSize)).copyTo(test);
	imwrite("64656_4.png", test*255.0);
	TrainLum(Rect(398, 102, TrainImgHRSize, TrainImgHRSize)).copyTo(test);
	imwrite("64656_5.png", test*255.0);
}
int main(int argc0, char **argv0) {
	int argc = argc0 - 1;
	char **argv = argv0 + 1;

	if (argc < 1) {
		usage();
	}
	if (!strcmp(argv[0], "Superes")) {
		assert2(argc == 4, "expected 4 args in Superes mode");
		string testname = argv[1];
		TRAINFILENAME = argv[2];
		string outputname = argv[3];
		LoadImg(testname, TRAINFILENAME, outputname);
	}
	else
	{
		assert2(argv[0] == "Superes", "expected only Superes mode");
	}

	//Test2();
	//Test1();
	GeneraTrainSet();
	SuperRe();
	waitKey(0);
	return 0;
}
