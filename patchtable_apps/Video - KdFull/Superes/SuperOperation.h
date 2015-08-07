#pragma once
#include "stdafx.h"
#include <flann/flann.hpp>
#include "patchtable.h"

#include <iostream>
#include <fstream>

using std::ifstream;
using std::stringstream;
using std::ofstream;
using std::cout;
using std::endl;

//#include <opencv2/imgproc/imgproc.hpp>
//using namespace cv;

/****
Input:
1.a pat of highresolution pic
2.lowresultion path
3. Upsample factor
4. SuperresResult filename
maybe use many pic for training data but the output result only one pic.
Suppose all the highresolution image is the same size.
**/
class SuperOperation
{
public:
	SuperOperation(string highresTxt, string lowresTxt, int Upsamplefactor, string outputfilename, int framePersInterval);
	~SuperOperation();
	void ParameterSetting();
	void ReadFilename();
	void LoadTrainImg(int startno, int number);
	void LoadTestImg(int testNo);
	void GenerateTrainSet();
	void Concatenat(cv::Mat lowpatch1, cv::Mat & random_vector, int keyrow);
	void Is_boundary(int &row, int &col, int rowpathsize, int colpathsize, int filter);
	void Filter(cv::Mat lowpatch1, cv::Mat & random_vector, int mode);
	void InWhichPic(int index, int lowCol, int &r, int &c, int& pic);
	float Average(cv::Mat& input);
	void ProcessImage(string ResultFile);
	
	int  In_whichPic(int &tempR, int rows);
	void PatchtableBuild();
	void PatchTableMethod();
	void KdTreeBuild();
	void KdTreeMethod();
	void BulidWeight();

	string im_highresTxt;
	string im_lowresTxt;
	string im_outputfilename;

	int im_DownSampleFactor = 0;
	int im_PCADIMS = 0;
	int im_picfortrainNumber=0;
	int TrainImgHRSize = 0;
	int TrainImgColsize = 0;
	int TrainImgRowsize = 0;
	int TrainImgLRSize = 0;
	int KeyRows = 0;
	int vectorSize = 0;
	int im_Kneighbour = 0;
	vector<cv::Mat> im_TrainLum;
	cv::Mat im_TestImg, im_TestLum, im_TestCR, im_TestCB;
	cv::Mat im_ResultImg;
	vector<string> TESTFILENAME;
	vector<string> TRAINFILENAME;
	//cv::Mat OriginalTrainFeatureData;
	cv::Mat randomKey;
	cv::Mat TrainFeatureData;
	int im_FramePersInterval;
	cv::PCA *im_pca;
	cv::Mat im_HighPatchSumWeight;
	/*PatchTableParams *p;
	PatchTable<>* im_table = NULL;
	cv::flann::Index im_flann_index;*/
	float *TrainFeatureArray;
	int TrainFeatureRows;
	int TrainFeatureCols;
	typedef flann::Index<flann::L2<float> > FlannIndexType;
	shared_ptr<FlannIndexType> flann_index;

#ifdef DrawResult
	int TrainMaxValue = 0;
	int TrainMinValue = 0;
#endif
};

