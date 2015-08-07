#include "SuperOperation.h"

SuperOperation::SuperOperation(string highresTxt, string lowresTxt, int Upsamplefactor, string outputfilename, int framePersInterval)
{
	im_highresTxt = highresTxt;
	im_lowresTxt = lowresTxt;
	im_DownSampleFactor = Upsamplefactor;
	im_outputfilename = outputfilename;
	TrainImgLRSize = 8;
	TrainImgHRSize = Upsamplefactor*TrainImgLRSize;
	KeyRows = 0;
	vectorSize = 256;
	im_FramePersInterval = framePersInterval;
	im_HighPatchSumWeight = cv::Mat::zeros(TrainImgHRSize, TrainImgHRSize, CV_32FC1);
	BulidWeight();
}


SuperOperation::~SuperOperation()
{
}
void SuperOperation::ParameterSetting()
{
	im_PCADIMS = 20;
	im_picfortrainNumber = 1;
	im_Kneighbour = 9;
}
void SuperOperation::BulidWeight()
{
	float *x = new  float[TrainImgHRSize];
	float *y = new  float[TrainImgHRSize];
	for (int i = 0; i < TrainImgHRSize; i++)
	{
		if (i < 9)
		{
			x[i] = 0.1 + 0.1*i;
		}
		else if (i > TrainImgHRSize - 10)
		{
			x[i] = 0.9 - 0.1*(9 - (TrainImgHRSize - i));
		}
		else
		{
			x[i] = 1.0;
		}
	}
	for (int i = 0; i < TrainImgHRSize; i++)
	{
		if (i < 9)
		{
			y[i] = 0.1 + 0.1*i;
		}
		else if (i > TrainImgHRSize - 10)
		{
			y[i] = 0.9 - 0.1*(9 - (TrainImgHRSize - i));
		}
		else
		{
			y[i] = 1.0;
		}
	}
	for (int i = 0; i < TrainImgHRSize; i++)
	{
		for (int j = 0; j < TrainImgHRSize; j++)
		{
			im_HighPatchSumWeight.at<float>(i, j) = x[i] * y[j];
		}
	}
	delete x;
	delete y;
}
void SuperOperation::Is_boundary(int &row, int &col, int rowpathsize, int colpathsize, int filter)
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
void SuperOperation::Concatenat(cv::Mat lowpatch1, cv::Mat & random_vector, int keyrow)
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
void SuperOperation::GenerateTrainSet()
{
	TrainImgColsize = im_TrainLum[0].cols;
	TrainImgRowsize = im_TrainLum[0].rows;
	KeyRows = (TrainImgColsize - TrainImgHRSize + 1)*(TrainImgRowsize - TrainImgHRSize + 1)*im_TrainLum.size();//highres patch number
	cout << "KeyRows:" << KeyRows << endl;
	cout << "train number:" << im_TrainLum.size() << endl;
	cv::RNG rng;
	int i = 0;
	int *temp = new int[KeyRows];
	randomKey.create(10000, vectorSize, CV_32FC1);
	for (int n = 0; n < KeyRows; n++)
	{
		temp[n] = 0;
	}
	TrainFeatureData = cv::Mat::zeros(KeyRows, im_PCADIMS, CV_32FC1);
	TrainFeatureRows = KeyRows;
	TrainFeatureCols = im_PCADIMS;
	TrainFeatureArray = new float[KeyRows*im_PCADIMS];
	for (int r = 0; r < KeyRows*im_PCADIMS; r++)
	{
		TrainFeatureArray[r] = 0.0;
	}
	std::ofstream LowResfile("time.txt", std::ofstream::out | std::ofstream::app);
	clock_t start, finish;
	double  duration = 0.0;
	start = clock();
	do{
		int col = rng.uniform(0, TrainImgColsize - TrainImgHRSize + 1);
		int row = rng.uniform(0, TrainImgRowsize - TrainImgHRSize + 1);
		int rngNumber = rng.uniform(0, im_TrainLum.size() - 1);
		if (temp[row*(TrainImgColsize - TrainImgHRSize + 1) + col + (TrainImgColsize - TrainImgHRSize + 1)*(TrainImgRowsize - TrainImgHRSize + 1)*rngNumber] == 0)
		{
			temp[row*(TrainImgColsize - TrainImgHRSize + 1) + col + (TrainImgColsize - TrainImgHRSize + 1)*(TrainImgRowsize - TrainImgHRSize + 1)*rngNumber] = 1;
			cv::Mat highpatch(TrainImgHRSize, TrainImgHRSize, CV_32FC1);
			cv::Mat lowpatch1 = cv::Mat::zeros(TrainImgLRSize, TrainImgLRSize, CV_32FC1);
			im_TrainLum[rngNumber](cv::Rect(col, row, TrainImgHRSize, TrainImgHRSize)).copyTo(highpatch);
			cv::resize(highpatch, lowpatch1, cv::Size(highpatch.cols / im_DownSampleFactor, highpatch.rows / im_DownSampleFactor), 0, 0, CV_INTER_CUBIC);
			Concatenat(lowpatch1, randomKey, i);
			i++;
		}
	} while (i < 10000);
	delete temp;

	im_pca = new  cv::PCA(randomKey, cv::Mat(), CV_PCA_DATA_AS_ROW, im_PCADIMS);
	
	int picRows = (TrainImgColsize - TrainImgHRSize + 1)*(TrainImgRowsize - TrainImgHRSize + 1);
	
#pragma omp parallel for
	for (int i = 0; i < im_TrainLum.size(); i++)
	{
		for (int r = 0; r < TrainImgRowsize - TrainImgHRSize + 1; r++)
		{
			for (int c = 0; c < TrainImgColsize - TrainImgHRSize + 1; c++)
			{
				cv::Mat input1(1, vectorSize, CV_32FC1);
				cv::Mat highpatch(TrainImgHRSize, TrainImgHRSize, CV_32FC1);
				cv::Mat lowpatch1 = cv::Mat::zeros(TrainImgLRSize, TrainImgLRSize, CV_32FC1);
				im_TrainLum[i](cv::Rect(c, r, TrainImgHRSize, TrainImgHRSize)).copyTo(highpatch);
				cv::resize(highpatch, lowpatch1, cv::Size(highpatch.cols / im_DownSampleFactor, highpatch.rows / im_DownSampleFactor), 0, 0, cv::INTER_AREA);
				Concatenat(lowpatch1, input1, 0);
				cv::Mat input2(1, im_PCADIMS, CV_32FC1);
				im_pca->project(input1, input2);
				for (int n = 0; n < input2.cols; n++)
				{
					TrainFeatureData.at<float>(r*(TrainImgColsize - TrainImgHRSize + 1) + c + i*picRows, n) = input2.at<float>(0, n);
					//TrainFeatureArray[(r*(TrainImgColsize - TrainImgHRSize + 1) + c + i*picRows)*im_PCADIMS + n] = input2.at<float>(0, n);
#ifdef DrawResult
					TrainMaxValue = (input1.at<float>(0, n) - TrainMaxValue)>0 ? input1.at<float>(0, n) : TrainMaxValue;
					TrainMinValue = (TrainMinValue - input1.at<float>(0, n))>0 ? input1.at<float>(0, n) : TrainMinValue;
#endif
				}
			}
		}
	}
	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	LowResfile << duration<< endl;
	LowResfile.close();
	for (int r = 0; r < KeyRows; r++)
	{
		for (int c = 0; c < im_PCADIMS; c++)
		{
			TrainFeatureArray[r*im_PCADIMS + c] =TrainFeatureData.at<float>(r, c);
		}
	}
}
void SuperOperation::Filter(cv::Mat lowpatch1, cv::Mat & random_vector, int mode)
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
void SuperOperation::LoadTrainImg(int startno, int number)
{
	cv::Mat tempTrain;
	vector< cv::Mat> LumYcYb;
	im_TrainLum.clear();
	for (int i = startno; i < startno+number; i++)
	{
		tempTrain = cv::imread(TRAINFILENAME[i], 1);
		tempTrain.convertTo(tempTrain, CV_32FC3);
		tempTrain = tempTrain / 255.0;
		cv::cvtColor(tempTrain, tempTrain, cv::COLOR_BGR2YCrCb);
		split(tempTrain, LumYcYb);
		im_TrainLum.push_back(LumYcYb[0].clone());
		LumYcYb.clear();
	}
}
void SuperOperation::LoadTestImg(int testNo)
{
	vector< cv::Mat> LumYcYb;
	im_TestImg = cv::imread(TESTFILENAME[testNo], 1);
	im_TestImg.convertTo(im_TestImg, CV_32FC3);
	im_TestImg = im_TestImg / 255.0;
	cv::cvtColor(im_TestImg, im_TestImg, cv::COLOR_BGR2YCrCb);
	split(im_TestImg, LumYcYb);
	im_TestLum = LumYcYb[0].clone();
	im_TestCR = LumYcYb[1].clone();
	im_TestCB = LumYcYb[2].clone();
	LumYcYb.clear();
}
void SuperOperation::InWhichPic(int index, int lowCol, int &r, int &c,int& pic)
{
	int n = 0;
	for (int i = 0; i < im_TrainLum.size();i++)
	{
		if ((index >= lowCol*(TrainImgRowsize - TrainImgHRSize + 1)*i) && index < lowCol*(TrainImgRowsize - TrainImgHRSize + 1)*(i + 1))
		{
		    pic = i;
			index = index - lowCol*(TrainImgRowsize - TrainImgHRSize + 1)*i;
			break;
		}
	}
	for (int x = 0; x < lowCol; x++)
	{
		if (index >= lowCol*x && index < lowCol*(x + 1))
		{
			r = x;
			c = index - x*lowCol;
			break;
		}

	}

}
void SuperOperation::ReadFilename()
{
	//read trainName
	int positionPotAfter = im_highresTxt.find_first_of(".");
	string testPotAfter = im_highresTxt.substr(0, positionPotAfter);
	string testForMat = im_highresTxt.substr(positionPotAfter + 1, im_highresTxt.size());
	if (testForMat == "txt")
	{
		ifstream in(im_highresTxt);
		char line[102];
		string nameTemp;
		stringstream ss2;
		while (in.getline(line, 102)){
			ss2 << line;
			ss2 >> nameTemp;
			TRAINFILENAME.push_back(nameTemp);
			ss2.clear();
		}
		in.close();
	}
	/*else
	{
		TESTFILENAME.push_back(testname);
	}*/
	//read testName
	positionPotAfter = im_lowresTxt.find_first_of(".");
	testPotAfter = im_lowresTxt.substr(0, positionPotAfter);
	testForMat = im_lowresTxt.substr(positionPotAfter + 1, im_lowresTxt.size());
	if (testForMat == "txt")
	{
		std::ifstream in(im_lowresTxt);
		char line[1024];
		string nameTemp;
		std::stringstream ss2;
		while (in.getline(line, 1024)){
			ss2 << line;
			ss2 >> nameTemp;
			TESTFILENAME.push_back(nameTemp);
			ss2.clear();
		}
		in.close();
	}
}
float SuperOperation::Average(cv::Mat& input)
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
int  SuperOperation::In_whichPic(int &tempR,int rows)
{
	int correctRow = 0,correctPic=0;
	for (int i = 0; i < im_TrainLum.size();i++)
	{
		if (tempR>=i*rows&&tempR < (i + 1)*rows)
		{
			tempR = tempR - i*rows;
			correctPic = i;
			break;
		}
	}
	return correctPic;
}
void SuperOperation::PatchtableBuild()
{
	//int trainh = (TrainImgRowsize - TrainImgHRSize + 1)*im_TrainLum.size();
	//int trainw = TrainImgColsize - TrainImgHRSize + 1;
	//Array<float> TrainVectorArray(trainh, trainw, im_PCADIMS);

	//int count = 0;
	//for (int r = 0; r < trainh; r++)
	//{
	//	for (int c = 0; c < trainw; c++)
	//	{
	//		for (int k = 0; k < im_PCADIMS; k++)
	//		{
	//			TrainVectorArray(r, c, k) = TrainFeatureData.at<float>(count, k);
	//		}
	//		count++;
	//	}
	//}
	//p = new PatchTableParams();
	//p->set_speed(0); // Needs to be above the other lines
	//p->coherence_spatial = 0.0;//1.2.3.4 or 0.5
	////p->calc_exact_dist = true;

	//p->grid_ndims = 10;
	//p->ndims = im_PCADIMS;
	//p->nchroma = 0;
	//p->is_descriptor = true;

	////p->limit = 4e6;    // Number of cells in table
	////p->do_rs = true;   // Whether to run random search
	////p->run_dt = true;  // Whether to fill in empty cells with distance transform
	////p->do_prop = true; // Whether to run propagation

	////allow contains the TrainFeatureData which is height*width rows and PCADIMS cols
	////ofstream LowResfile("time.txt", std::ofstream::out | std::ofstream::app);
	////clock_t start, finish;
	////double  duration;
	////start = clock();
	//im_table = new PatchTable<>(p, TrainVectorArray);
	////finish = clock();
	////duration = (double)(finish - start) / CLOCKS_PER_SEC;
	////LowResfile << "PatchTable build" << duration << ".seconds" << endl;
	//TrainFeatureData.release();
}
void SuperOperation::KdTreeBuild()
{
	std::ofstream LowResfile("time.txt", std::ofstream::out | std::ofstream::app);
	clock_t start, finish;
	double  duration;
	start = clock();
	//////////////////////////////FOR OPENCV FLANNKDTreeIndexParams：4 before for excel
	// im_flann_index.build(TrainFeatureData,
    //  flann::KDTreeIndexParams(1),
	//	cvflann::FLANN_DIST_L2
	//	);
	flann::Matrix<float> train_data(TrainFeatureArray, TrainFeatureRows, TrainFeatureCols);
	flann_index = make_shared<FlannIndexType>(train_data, flann::KDTreeSingleIndexParams(10));
	flann_index->buildIndex();
	/*float speedup;
	p = DEFAULT_FLANN_PARAMETERS;
	p.algorithm = FLANN_INDEX_KDTREE;
	p.trees = 4;
	p.log_level = FLANN_LOG_INFO;
	p.checks = 8;
	p.eps = 50.0;
	index_id = flann_build_index(TrainFeatureData, TrainFeatureRows, TrainFeatureCols, &speedup, &p);*/
	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	LowResfile << duration  << std::endl;
	LowResfile.close();
	//delete TrainFeatureData;
	//TrainFeatureData.release();
}
void SuperOperation::PatchTableMethod()
{
//		int lowRows = im_TestImg.rows;
//		int lowCols = im_TestImg.cols;
//		int tcount = (lowRows - TrainImgLRSize + 1)*(lowCols - TrainImgLRSize + 1);
//		cv::Mat filter1mat = cv::Mat::zeros(lowRows, lowCols, CV_32FC1);
//		cv::Mat filter2mat = cv::Mat::zeros(lowRows, lowCols, CV_32FC1);
//		cv::Mat filter3mat = cv::Mat::zeros(lowRows, lowCols, CV_32FC1);
//		cv::Mat filter4mat = cv::Mat::zeros(lowRows, lowCols, CV_32FC1);
//		Filter(im_TestLum, filter1mat, 1);
//		Filter(im_TestLum, filter2mat, 2);
//		Filter(im_TestLum, filter3mat, 3);
//		Filter(im_TestLum, filter4mat, 4);
//		im_ResultImg = cv::Mat::zeros(lowRows*im_DownSampleFactor, lowCols*im_DownSampleFactor, CV_32FC3);
//		cv::Mat finalresult = cv::Mat::zeros(lowRows*im_DownSampleFactor, lowCols*im_DownSampleFactor, CV_32FC3);
//		cv::Mat SplitTest = cv::Mat::zeros((lowRows - TrainImgLRSize + 1)*TrainImgLRSize, (lowCols - TrainImgLRSize + 1)*TrainImgLRSize, CV_32FC1);
//		cv::Mat testemp(TrainImgLRSize, TrainImgLRSize, CV_32FC1);
//		int coutcount = 0;
//		cv::Mat TestFeatureData = cv::Mat::zeros(tcount, im_PCADIMS, CV_32FC1);
//		cv::Mat TestFeatureDataOriginal = cv::Mat::zeros(tcount, 256, CV_32FC1);
//
//		cv::Mat GaussHighPatchWeight = cv::Mat::zeros(lowRows - TrainImgLRSize + 1, lowCols - TrainImgLRSize + 1, CV_32FC1);
//		for (int r = 0; r < lowRows - TrainImgLRSize + 1; r++)
//		{
//			for (int c = 0; c < lowCols - TrainImgLRSize + 1; c++)
//			{
//				//r = 6, c = 11;
//				cv::Mat input1(1, vectorSize, CV_32FC1);
//				cv::Mat output1(1, im_PCADIMS, CV_32FC1);
//				cv::Mat testMatch = cv::Mat::ones(TrainImgHRSize + 1, (8 + 1 + 8 + 1 + 8 + 1 + 8 + 72 + 1), CV_32FC1);
//				im_TestLum(cv::Rect(c, r, TrainImgLRSize, TrainImgLRSize)).copyTo(testemp);
//				testemp.copyTo(SplitTest(cv::Rect(c*TrainImgLRSize, r*TrainImgLRSize, TrainImgLRSize, TrainImgLRSize)));
//				filter1mat(cv::Rect(c, r, TrainImgLRSize, TrainImgLRSize)).copyTo(testemp);
//				int n = 0;
//				for (int i = 0; i < TrainImgLRSize; i++)
//				{
//					for (int j = 0; j < TrainImgLRSize; j++)
//					{
//						//input1.at<float>(0, n) = (testemp.at<float>(i, j) - MinValue) / (MaxValue - MinValue);
//						input1.at<float>(0, n) = testemp.at<float>(i, j);
//						n++;
//					}
//				}
//				filter2mat(cv::Rect(c, r, TrainImgLRSize, TrainImgLRSize)).copyTo(testemp);
//				for (int i = 0; i < TrainImgLRSize; i++)
//				{
//					for (int j = 0; j < TrainImgLRSize; j++)
//					{
//						//input1.at<float>(0, n) = (testemp.at<float>(i, j) - MinValue) / (MaxValue - MinValue);
//						input1.at<float>(0, n) = testemp.at<float>(i, j);
//						n++;
//					}
//				}
//				filter3mat(cv::Rect(c, r, TrainImgLRSize, TrainImgLRSize)).copyTo(testemp);
//				for (int i = 0; i < TrainImgLRSize; i++)
//				{
//					for (int j = 0; j < TrainImgLRSize; j++)
//					{
//						//input1.at<float>(0, n) = (testemp.at<float>(i, j) - MinValue) / (MaxValue - MinValue);
//						input1.at<float>(0, n) = testemp.at<float>(i, j);
//						n++;
//					}
//				}
//				filter4mat(cv::Rect(c, r, TrainImgLRSize, TrainImgLRSize)).copyTo(testemp);
//				for (int i = 0; i < TrainImgLRSize; i++)
//				{
//					for (int j = 0; j < TrainImgLRSize; j++)
//					{
//						//input1.at<float>(0, n) = (testemp.at<float>(i, j) - MinValue) / (MaxValue - MinValue);
//						input1.at<float>(0, n) = testemp.at<float>(i, j);
//						n++;
//					}
//				}
//				im_pca->project(input1, output1);
//				for (int q = 0; q < output1.cols; q++)
//				{
//					TestFeatureData.at<float>(coutcount, q) = output1.at<float>(0, q);
//				}
//				for (int p = 0; p < input1.cols; p++)
//				{
//					TestFeatureDataOriginal.at<float>(coutcount, p) = input1.at<float>(0, p);
//				}
//				coutcount++;
//			}
//		}
//		
//		int testh = lowRows - TrainImgLRSize + 1;
//		int testw = lowCols - TrainImgLRSize + 1;
//		Array<float> TestVectorArray(testh, testw, im_PCADIMS);
//		int count = 0;
//		for (int r = 0; r < testh; r++)
//		{
//			for (int c = 0; c < testw; c++)
//			{
//				for (int k = 0; k < im_PCADIMS; k++)
//				{
//					TestVectorArray(r, c, k) = TestFeatureData.at<float>(count, k);
//				}
//				count++;
//			}
//		}
//		Array<double> ann;
//		double latest_time = 1e100;
//		for (int iter = 0; iter < 1; iter++){
//			//start = clock();
//			latest_time = im_table->lookup(TestVectorArray, ann);
//			//finish = clock();
//		}
//		//duration = (double)(finish - start) / CLOCKS_PER_SEC;
//		//LowResfile << "look up method" << duration << ".seconds" << endl;
//		//LowResfile.close();
//
//		//ReaDate(KeyRows, vectorSize, testset, "testset.txt");
//		//ReaDate(tcount, nn, result, "result.txt");
//		/********************Compute Weight***/
//		
//		cv::Mat weightResult = cv::Mat::zeros(lowRows*im_DownSampleFactor, lowCols*im_DownSampleFactor, CV_32FC1);
//		cv::Mat Result = cv::Mat::zeros(lowRows*im_DownSampleFactor, lowCols*im_DownSampleFactor, CV_32FC1);
//		int r = 0, c = 0, tempR = 0, tempC = 0; int index[9];
//		//cv::Mat resultFindHighPatch = cv::Mat::zeros(TrainImgHRSize * 50, TrainImgHRSize*(nn), CV_32FC1);
//		//cv::Mat resultFindLowPatch = cv::Mat::zeros(TrainImgLRSize * 50, TrainImgLRSize*(nn + 1), CV_32FC1);
//		cv::Mat tempHighPatchSum = cv::Mat::zeros(TrainImgHRSize, TrainImgHRSize, CV_32FC1);
//
//		int ttcount = 0, correctPic=0;
//		for (int r = 0; r < (lowRows - TrainImgLRSize + 1); r++)
//		{
//			for (int c = 0; c < (lowCols - TrainImgLRSize + 1); c++)
//			{
//				cv::Mat tempHighPatch = cv::Mat::zeros(TrainImgHRSize, TrainImgHRSize, CV_32FC1);
//				cv::Mat tempLowPatch = cv::Mat::zeros(TrainImgLRSize, TrainImgLRSize, CV_32FC1);
//				cv::Mat input1(1, vectorSize, CV_32FC1);
//				cv::Mat highpatch(TrainImgHRSize, TrainImgHRSize, CV_32FC1);
//				cv::Mat lowpatch1 = cv::Mat::zeros(TrainImgLRSize, TrainImgLRSize, CV_32FC1);
//				im_TestLum(Range(r, r + TrainImgLRSize), Range(c, c + TrainImgLRSize)).copyTo(tempLowPatch);
//				for (int i = 0; i < 1; i++)
//				{
//					tempC = ann(r, c, NNF_X);
//					tempR = ann(r, c, NNF_Y);
//					
//					/*if (tempR > TrainImgRowsize)
//					{
//						cout << "error: trainLum's cv::Size is more than one" << endl;
//					}*/
//					correctPic = In_whichPic(tempR, (TrainImgRowsize - TrainImgHRSize + 1));
//					float sum = 0.0, tempsum = 0.0;
//					im_TrainLum[correctPic](cv::Rect(tempC, tempR, TrainImgHRSize, TrainImgHRSize)).copyTo(highpatch);
//					cv::resize(highpatch, lowpatch1, cv::Size(highpatch.cols / im_DownSampleFactor, highpatch.rows / im_DownSampleFactor), 0, 0, INTER_AREA);
//					Concatenat(lowpatch1, input1, 0);
//					for (int x = 0; x<256; x++)
//					{
//						tempsum = abs(input1.at<float>(0, x) - TestFeatureDataOriginal.at<float>(r*(lowCols - TrainImgLRSize + 1) + c, x));
//						sum = tempsum*tempsum + sum;
//					}
//					sum = exp(-sum / (2 * 40));
//					GaussHighPatchWeight.at<float>(r, c) = sum;
//#ifdef DrawResult
//					cv::Mat input1(1, vectorSize, CV_32FC1);
//					filter1mat(cv::Rect(c, r, TrainImgLRSize, TrainImgLRSize)).copyTo(testemp);
//					int n = 0;
//					for (int i = 0; i < TrainImgLRSize; i++)
//					{
//						for (int j = 0; j < TrainImgLRSize; j++)
//						{
//							//input1.at<float>(0, n) = (testemp.at<float>(i, j) - MinValue) / (MaxValue - MinValue);
//							input1.at<float>(0, n) = testemp.at<float>(i, j);
//							n++;
//						}
//					}
//					filter2mat(cv::Rect(c, r, TrainImgLRSize, TrainImgLRSize)).copyTo(testemp);
//					for (int i = 0; i < TrainImgLRSize; i++)
//					{
//						for (int j = 0; j < TrainImgLRSize; j++)
//						{
//							//input1.at<float>(0, n) = (testemp.at<float>(i, j) - MinValue) / (MaxValue - MinValue);
//							input1.at<float>(0, n) = testemp.at<float>(i, j);
//							n++;
//						}
//					}
//					filter3mat(cv::Rect(c, r, TrainImgLRSize, TrainImgLRSize)).copyTo(testemp);
//					for (int i = 0; i < TrainImgLRSize; i++)
//					{
//						for (int j = 0; j < TrainImgLRSize; j++)
//						{
//							//input1.at<float>(0, n) = (testemp.at<float>(i, j) - MinValue) / (MaxValue - MinValue);
//							input1.at<float>(0, n) = testemp.at<float>(i, j);
//							n++;
//						}
//					}
//					filter4mat(cv::Rect(c, r, TrainImgLRSize, TrainImgLRSize)).copyTo(testemp);
//					for (int i = 0; i < TrainImgLRSize; i++)
//					{
//						for (int j = 0; j < TrainImgLRSize; j++)
//						{
//							//input1.at<float>(0, n) = (testemp.at<float>(i, j) - MinValue) / (MaxValue - MinValue);
//							input1.at<float>(0, n) = testemp.at<float>(i, j);
//							n++;
//						}
//					}
//					cv::Mat testMatch = cv::Mat::ones(TrainImgHRSize + 1, (8 + 1 + 8 + 1 + 8 + 1 + 8 + 72 + 1), CV_32FC1);
//					TestLum(cv::Rect(c, r, TrainImgLRSize, TrainImgLRSize)).copyTo(testemp);
//					testemp.copyTo(SplitTest(cv::Rect(c*TrainImgLRSize, r*TrainImgLRSize, TrainImgLRSize, TrainImgLRSize)));
//					int testr = 0, testc = 0; //float distance = 0;
//					for (int j = 0; j < vectorSize; j++) {
//						testMatch.at<float>(testr, testc) = (input1.at<float>(0, j) - MinValue) / (MaxValue - MinValue);
//						testMatch.at<float>(testr, testc + 9) = (OriginalTrainFeatureData.at<float>(tempR*(TrainImgColsize - TrainImgHRSize + 1) + tempC, j) - TrainMinValue) / (TrainMaxValue - TrainMinValue);
//						//cout << testMatch.at<float>(testr, testc) << " " << testMatch.at<float>(testr, testc + 9) << " " << testMatch.at<float>(testr, testc) - testMatch.at<float>(testr, testc + 9)<<endl;
//						//distance = abs(testMatch.at<float>(testr, testc) - testMatch.at<float>(testr, testc + 9))+distance;
//						testc++;
//						if (testc == 8)
//						{
//							testc = 0; testr++;
//						}
//					}
//					//cout << "the distance of the original is" << 13.4538 << ";   the distance of the our result is" << distance;
//					for (int tesr = 0; tesr < TrainImgLRSize; tesr++)
//					{
//						for (int tesc = 0; tesc < TrainImgLRSize; tesc++)
//						{
//							testMatch.at<float>(tesr, 18 + tesc) = TestLum.at<float>(r + tesr, c + tesc);
//						}
//					}
//					cv::Mat testemphr(TrainImgHRSize, TrainImgHRSize, CV_32FC1);
//					for (int r = 0; r < TrainImgHRSize; r++)
//					{
//						for (int c = 0; c < TrainImgHRSize; c++)
//						{
//							testMatch.at<float>(0 + r, 36 + c) = TrainLum.at<float>(tempR + r, tempC + c);
//							testemphr.at<float>(r, c) = TrainLum.at<float>(tempR + r, tempC + c);
//						}
//					}
//					resize(testemphr, testemp, cv::Size(TrainImgLRSize, TrainImgLRSize), 0, 0, INTER_CUBIC);
//					for (int r = 0; r < TrainImgLRSize; r++)
//					{
//						for (int c = 0; c < TrainImgLRSize; c++)
//						{
//							testMatch.at<float>(0 + r, 27 + c) = testemp.at<float>(r, c);
//						}
//					}
//					char file_name[15];
//					sprintf(file_name, "TestFirstMatch%d.png", tempR*(TrainImgColsize - TrainImgHRSize + 1) + tempC);
//					cv::imwrite(file_name, testMatch*255.0);
//#endif
//					float average_highpatch = Average(highpatch);
//					float average_lowpatch = Average(tempLowPatch);
//					//cout << average_highpatch << endl;
//					//cout << average_lowpatch << endl;
//					//resize(tempHighPatch, tempLowPatch, cv::Size(tempHighPatch.cols / DownSampleFactor, tempHighPatch.rows / DownSampleFactor), 0, 0, CV_INTER_CUBIC);
//					//tempHighPatchSum = exp(param) * tempHighPatch + tempHighPatchSum;
//					tempHighPatchSum = highpatch - average_highpatch + average_lowpatch;
//
//				}
//				//imwrite("TestFinal_11.png", highpatch*255.0);
//				//char file_name[15];
//				//sprintf(file_name, "TestFinal%d.png", ttcount);
//				//imwrite(file_name, tempHighPatchSum*255.0);
//				ttcount++;
//				int startIndexR = r * im_DownSampleFactor; int startIndexC = c * im_DownSampleFactor;
//				for (int i = 0; i < TrainImgHRSize; i++)
//				{
//					for (int j = 0; j < TrainImgHRSize; j++)
//					{
//						Result.at<float>(i + startIndexR, j + startIndexC) = tempHighPatchSum.at<float>(i, j)*im_HighPatchSumWeight.at<float>(i, j)*GaussHighPatchWeight.at<float>(r, c) + Result.at<float>(i + startIndexR, j + startIndexC);
//						//Result.at<float>(i + startIndexR, j + startIndexC) = tempHighPatchSum.at<float>(i, j) + Result.at<float>(i + startIndexR, j + startIndexC);
//						//weightResult.at<float>(i + startIndexR, j + startIndexC) = im_HighPatchSumWeight.at<float>(i, j) + weightResult.at<float>(i + startIndexR, j + startIndexC);
//						weightResult.at<float>(i + startIndexR, j + startIndexC) = im_HighPatchSumWeight.at<float>(i, j)*GaussHighPatchWeight.at<float>(r, c) + weightResult.at<float>(i + startIndexR, j + startIndexC);
//					}
//				}
//			}
//		}
//		for (int i = 0; i < Result.rows; i++)
//		{
//			for (int j = 0; j < Result.cols; j++)
//			{
//				Result.at<float>(i, j) = Result.at<float>(i, j) / weightResult.at<float>(i, j);
//			}
//		}
//		//char file_name[50];
//		//sprintf(file_name, "addOriginalResult%d.png", im_PCADIMS);
//		//imwrite(TEST + "addOriginalResult" + PCADIMS + ".png", Result*255.0);
//		//imwrite(file_name, Result*255.0);
//		//im_TestLum = Result;
//		//imwrite(file_name, Result*255.0);
//		cv::Mat CR = Result.clone();
//		cv::Mat CB = Result.clone();
//		cv::Mat CRt = cv::Mat::zeros(CR.rows, CR.cols, CV_32FC1);
//		cv::Mat CBt = cv::Mat::zeros(CR.rows, CR.cols, CV_32FC1);
//		cv::resize(im_TestCR, CR, cv::Size(im_TestCB.cols * im_DownSampleFactor, im_TestCB.rows * im_DownSampleFactor), 0, 0, INTER_CUBIC);
//		cv::bilateralFilter(CR, CRt, -1.0, 0.1, im_DownSampleFactor + 2);
//		cv::resize(im_TestCB, CB, cv::Size(im_TestCB.cols * im_DownSampleFactor, im_TestCB.rows * im_DownSampleFactor), 0, 0, INTER_CUBIC);
//		cv::bilateralFilter(CB, CBt, -1.0, 0.1, im_DownSampleFactor + 2);
//		vector< cv::Mat> LumYcYb;
//		LumYcYb.push_back(Result.clone());
//		LumYcYb.push_back(CRt.clone());
//		LumYcYb.push_back(CBt.clone());
//		cv::merge(LumYcYb, im_ResultImg);
//		cv::cvtColor(im_ResultImg, im_ResultImg, COLOR_YCrCb2BGR);
//		im_ResultImg = im_ResultImg * 255.0;
//		LumYcYb.clear();
}
void SuperOperation::KdTreeMethod()
{
	int index = 0, resultIndex = 0;
	
	//vector<int> indices(im_Kneighbour);
	//vector<float> dists(im_Kneighbour);

#ifdef DrawResult
	float MaxValue = -10000.0, MinValue = 100000.0;
#endif
	
		int lowRows = im_TestLum.rows;
		int lowCols = im_TestLum.cols;
		int tcount = (lowRows - TrainImgLRSize + 1)*(lowCols - TrainImgLRSize + 1);
		int **result = new int*[tcount];
		for (int n = 0; n < tcount; n++)
		{
			result[n] = new int[im_Kneighbour];
		}
		for (int n = 0; n < tcount; n++)
		{
			for (int y = 0; y < im_Kneighbour; y++)
			{
				result[n][y] = 0;
			}
		}
		cv::Mat filter1mat = cv::Mat::zeros(lowRows, lowCols, CV_32FC1);
		cv::Mat filter2mat = cv::Mat::zeros(lowRows, lowCols, CV_32FC1);
		cv::Mat filter3mat = cv::Mat::zeros(lowRows, lowCols, CV_32FC1);
		cv::Mat filter4mat = cv::Mat::zeros(lowRows, lowCols, CV_32FC1);
		Filter(im_TestLum, filter1mat, 1);
		Filter(im_TestLum, filter2mat, 2);
		Filter(im_TestLum, filter3mat, 3);
		Filter(im_TestLum, filter4mat, 4);
		im_ResultImg = cv::Mat::zeros(lowRows*im_DownSampleFactor, lowCols*im_DownSampleFactor, CV_32FC3);
		cv::Mat finalresult = cv::Mat::zeros(lowRows*im_DownSampleFactor, lowCols*im_DownSampleFactor, CV_32FC3);
		//cv::Mat SplitTest = cv::Mat::zeros((lowRows - TrainImgLRSize + 1)*TrainImgLRSize, (lowCols - TrainImgLRSize + 1)*TrainImgLRSize, CV_32FC1);
		cv::Mat testemp(TrainImgLRSize, TrainImgLRSize, CV_32FC1);
		int coutcount = 0;
#ifdef Patchtable
		cv::Mat TestFeatureData = cv::Mat::zeros(tcount, PCADIMS, CV_32FC1);
#endif
#ifdef DrawResult
		for (int r = 0; r < lowRows - TrainImgLRSize + 1; r++)
		{
			for (int c = 0; c < lowCols - TrainImgLRSize + 1; c++)
			{
				cv::Mat input1(1, vectorSize, CV_32FC1);
				cv::Mat output1(1, 20, CV_32FC1);
				cv::Mat testMatch = cv::Mat::ones(TrainImgHRSize + 1, (8 + 1 + 8 + 1 + 8 + 1 + 8 + 72 + 1), CV_32FC1);
				TestLum(cv::Rect(c, r, TrainImgLRSize, TrainImgLRSize)).copyTo(testemp);
				testemp.copyTo(SplitTest(cv::Rect(c*TrainImgLRSize, r*TrainImgLRSize, TrainImgLRSize, TrainImgLRSize)));
				filter1mat(cv::Rect(c, r, TrainImgLRSize, TrainImgLRSize)).copyTo(testemp);
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
				filter2mat(cv::Rect(c, r, TrainImgLRSize, TrainImgLRSize)).copyTo(testemp);
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
				filter3mat(cv::Rect(c, r, TrainImgLRSize, TrainImgLRSize)).copyTo(testemp);
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
				filter4mat(cv::Rect(c, r, TrainImgLRSize, TrainImgLRSize)).copyTo(testemp);
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
		std::ofstream LowResfile("time.txt", std::ofstream::out | std::ofstream::app);
		clock_t start, finish;
		double duration = 0.0;
		cv::Mat GaussHighPatchWeight = cv::Mat::zeros(lowRows - TrainImgLRSize + 1, lowCols - TrainImgLRSize + 1, CV_32FC1);
		cv::Mat TestFeatureDataOriginal = cv::Mat::zeros(tcount, 256, CV_32FC1);
		cv::Mat highpatch(TrainImgHRSize, TrainImgHRSize, CV_32FC1);
		cv::Mat lowpatch1 = cv::Mat::zeros(TrainImgLRSize, TrainImgLRSize, CV_32FC1);
		float *TestFeatureDataFloat = new float[im_PCADIMS];
		for (int r = 0; r < im_PCADIMS; r++)
		{
			TestFeatureDataFloat[r] = 0.0;
		}
		for (int r = 0; r < lowRows - TrainImgLRSize + 1; r++)
		{
			for (int c = 0; c < lowCols - TrainImgLRSize + 1; c++)
			{
				//r = 6, c = 11;
				cv::Mat input1(1, vectorSize, CV_32FC1);
				cv::Mat output1(1, im_PCADIMS, CV_32FC1);
				cv::Mat testMatch = cv::Mat::ones(TrainImgHRSize + 1, (8 + 1 + 8 + 1 + 8 + 1 + 8 + 72 + 1), CV_32FC1);
				//im_TestImg(cv::Rect(c, r, TrainImgLRSize, TrainImgLRSize)).copyTo(testemp);
				//testemp.copyTo(SplitTest(cv::Rect(c*TrainImgLRSize, r*TrainImgLRSize, TrainImgLRSize, TrainImgLRSize)));
				filter1mat(cv::Rect(c, r, TrainImgLRSize, TrainImgLRSize)).copyTo(testemp);
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
				filter2mat(cv::Rect(c, r, TrainImgLRSize, TrainImgLRSize)).copyTo(testemp);
				for (int i = 0; i < TrainImgLRSize; i++)
				{
					for (int j = 0; j < TrainImgLRSize; j++)
					{
						//input1.at<float>(0, n) = (testemp.at<float>(i, j) - MinValue) / (MaxValue - MinValue);
						input1.at<float>(0, n) = testemp.at<float>(i, j);
						n++;
					}
				}
				filter3mat(cv::Rect(c, r, TrainImgLRSize, TrainImgLRSize)).copyTo(testemp);
				for (int i = 0; i < TrainImgLRSize; i++)
				{
					for (int j = 0; j < TrainImgLRSize; j++)
					{
						//input1.at<float>(0, n) = (testemp.at<float>(i, j) - MinValue) / (MaxValue - MinValue);
						input1.at<float>(0, n) = testemp.at<float>(i, j);
						n++;
					}
				}
				filter4mat(cv::Rect(c, r, TrainImgLRSize, TrainImgLRSize)).copyTo(testemp);
				for (int i = 0; i < TrainImgLRSize; i++)
				{
					for (int j = 0; j < TrainImgLRSize; j++)
					{
						//input1.at<float>(0, n) = (testemp.at<float>(i, j) - MinValue) / (MaxValue - MinValue);
						input1.at<float>(0, n) = testemp.at<float>(i, j);
						n++;
					}
				}
				im_pca->project(input1, output1);
				
				for (int p = 0; p < input1.cols; p++)
				{
					TestFeatureDataOriginal.at<float>(coutcount, p) = input1.at<float>(0, p);
				}
				for (int q = 0; q < output1.cols; q++)
				{
					TestFeatureDataFloat[q] = output1.at<float>(0, q);
				}
#ifdef Patchtable
				for (int q = 0; q < output1.cols; q++)
				{
					TestFeatureData.at<float>(coutcount, q) = output1.at<float>(0, q);
				}
#else
				start = clock();
				//im_flann_index.knnSearch(output1, indices, dists, im_Kneighbour, flann::SearchParams(6));//SearchParams 256
				/*int* resultflann = (int*)malloc(1*im_Kneighbour*sizeof(int));
				float* dists = (float*)malloc(1*im_Kneighbour*sizeof(float));*/
				//index_id->knnSearch(query_matrix, matched_indices_mat, matched_dists_mat, knn, search_params);
				//flann_find_nearest_neighbors_index(index_id, TestFeatureDataFloat, 1, resultflann, dists, im_Kneighbour, &p);
				flann::Matrix<int> resultflann(new int[1 * im_Kneighbour], 1, im_Kneighbour);
				flann::Matrix<float> dists(new float[1 * im_Kneighbour], 1, im_Kneighbour);

				flann::SearchParams search_params(1,1600);
				//search_params.eps = 50.0;
				flann::Matrix<float> test_data(TestFeatureDataFloat, 1, im_PCADIMS);
				flann_index->knnSearch(test_data, resultflann, dists, im_Kneighbour, search_params);
				//flann_find_nearest_neighbors(TrainFeatureData, TrainFeatureRows, TrainFeatureCols, TestFeatureDataFloat, 1, resultflann, dists, im_Kneighbour,&p);
				finish = clock();
				duration = (double)(finish - start) / CLOCKS_PER_SEC + duration;
				//vector<int>::iterator indexn = indices.begin(); 
				int x2 = 0; int matchr = 0, matchc = 0;
				//vector<float>::iterator indexndis = dists.begin();
				float disMin = 10000.0;
				for (int i = 0; i < im_Kneighbour;i++)//for (; indexn != indices.end(); indexn++)
				{
					int n= resultflann[0][i];//int n = *indexn;
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
							testMatch.at<float>(tesr, 18 + tesc) = TestLum.at<float>(r + tesr, c + tesc);
						}
					}
					cv::Mat testemphr(TrainImgHRSize, TrainImgHRSize, CV_32FC1);
					for (int r = 0; r < TrainImgHRSize; r++)
					{
						for (int c = 0; c < TrainImgHRSize; c++)
						{
							testMatch.at<float>(0 + r, 36 + c) = TrainLum.at<float>(matchr + r, matchc + c);
							testemphr.at<float>(r, c) = TrainLum.at<float>(matchr + r, matchc + c);
						}
					}
					resize(testemphr, testemp, cv::Size(TrainImgLRSize, TrainImgLRSize), 0, 0, INTER_CUBIC);
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
					if (dists[0][i]<disMin)//if (*indexndis < disMin)
					{
						disMin = dists[0][i];//disMin = *indexndis;
						result[coutcount][x2] = result[coutcount][0]; x2++;
						result[coutcount][0] = n;
					}
					else
					{
						result[coutcount][x2] = n; x2++;
					}
					//indexndis++;
				}
#endif
				coutcount++;
			}
		}
		delete TestFeatureDataFloat;
		LowResfile << duration  << std::endl;
		LowResfile.close();
#ifndef Patchtable
		//cout << "search time:" << duration << ".seconds";
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
					TrainVectorArray(r, c, k) = TrainFeatureData.at<float>(count, k);
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
		/*ReaDate(KeyRows, vectorSize, testset, "testset.txt");
		/*ReaDate(tcount, nn, result, "result.txt");
		/********************Compute Weight**
		/*float *x = new  float[TrainImgHRSize];
		float *y = new  float[TrainImgHRSize];
		for (int i = 0; i < TrainImgHRSize; i++)
		{
			if (i < 9)
			{
				x[i] = 0.1 + 0.1*i;
			}
			else if (i > TrainImgHRSize - 10)
			{
				x[i] = 0.9 - 0.1*(9 - (TrainImgHRSize - i));
			}
			else
			{
				x[i] = 1.0;
			}
		}
		for (int i = 0; i < TrainImgHRSize; i++)
		{
			if (i < 9)
			{
				y[i] = 0.1 + 0.1*i;
			}
			else if (i > TrainImgHRSize - 10)
			{
				y[i] = 0.9 - 0.1*(9 - (TrainImgHRSize - i));
			}
			else
			{
				y[i] = 1.0;
			}
		}*/
		/*cv::Mat tempHighPatchWeight = cv::Mat::zeros(TrainImgHRSize, TrainImgHRSize, CV_32FC1);
		/*float **result = new float*[TrainImgHRSize];
		for (int n = 0; n < TrainImgHRSize; n++)
		{
			result[n] = new float[TrainImgHRSize];
		}
		for (int i = 0; i < TrainImgHRSize; i++)
		{
			for (int j = 0; j < TrainImgHRSize; j++)
			{
				tempHighPatchWeight.at<float>(i, j) = x[i] * y[j];
				result[i][j] = tempHighPatchWeight.at<float>(i, j);
			}
		}
		write_results("weight.txt", result, TrainImgHRSize, TrainImgHRSize);*/
		cv::Mat weightResult = cv::Mat::zeros(lowRows*im_DownSampleFactor, lowCols*im_DownSampleFactor, CV_32FC1);
		cv::Mat Result = cv::Mat::zeros(lowRows*im_DownSampleFactor, lowCols*im_DownSampleFactor, CV_32FC1);
		int r = 0, c = 0, tempR = 0, tempC = 0; float sumweight = 0.0; 
		cv::Mat resultFindHighPatch = cv::Mat::zeros(TrainImgHRSize * 50, TrainImgHRSize*(im_Kneighbour), CV_32FC1);
		cv::Mat resultFindLowPatch = cv::Mat::zeros(TrainImgLRSize * 50, TrainImgLRSize*(im_Kneighbour + 1), CV_32FC1);
		
		int ttcount = 0;
		for (int r = 0; r < (lowRows - TrainImgLRSize + 1); r++)
		{
			for (int c = 0; c < (lowCols - TrainImgLRSize + 1); c++)
			{
				//r = 6, c = 11;
				cv::Mat tempHighPatchSum = cv::Mat::zeros(TrainImgHRSize, TrainImgHRSize, CV_32FC1);
				cv::Mat tempHighPatch = cv::Mat::zeros(TrainImgHRSize, TrainImgHRSize, CV_32FC1);
				cv::Mat tempLowPatch = cv::Mat::zeros(TrainImgLRSize, TrainImgLRSize, CV_32FC1);
				im_TestLum(cv::Range(r, r + TrainImgLRSize), cv::Range(c, c + TrainImgLRSize)).copyTo(tempLowPatch);
				cv::Mat input1(1, vectorSize, CV_32FC1);
				cv::Mat lowpatch1 = cv::Mat::zeros(TrainImgLRSize, TrainImgLRSize, CV_32FC1);
				for (int i = 0; i < 1; i++)
				{
#ifndef	Patchtable
					index = result[ttcount][i];
					int pic=0;
					InWhichPic(index, (TrainImgColsize - TrainImgHRSize + 1), tempR, tempC,pic);
#else
					tempC = an
						n(r, c, NNF_X);
					tempR = ann(r, c, NNF_Y);
#ifdef DrawResult
					cv::Mat input1(1, vectorSize, CV_32FC1);
					filter1mat(cv::Rect(c, r, TrainImgLRSize, TrainImgLRSize)).copyTo(testemp);
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
					filter2mat(cv::Rect(c, r, TrainImgLRSize, TrainImgLRSize)).copyTo(testemp);
					for (int i = 0; i < TrainImgLRSize; i++)
					{
						for (int j = 0; j < TrainImgLRSize; j++)
						{
							//input1.at<float>(0, n) = (testemp.at<float>(i, j) - MinValue) / (MaxValue - MinValue);
							input1.at<float>(0, n) = testemp.at<float>(i, j);
							n++;
						}
					}
					filter3mat(cv::Rect(c, r, TrainImgLRSize, TrainImgLRSize)).copyTo(testemp);
					for (int i = 0; i < TrainImgLRSize; i++)
					{
						for (int j = 0; j < TrainImgLRSize; j++)
						{
							//input1.at<float>(0, n) = (testemp.at<float>(i, j) - MinValue) / (MaxValue - MinValue);
							input1.at<float>(0, n) = testemp.at<float>(i, j);
							n++;
						}
					}
					filter4mat(cv::Rect(c, r, TrainImgLRSize, TrainImgLRSize)).copyTo(testemp);
					for (int i = 0; i < TrainImgLRSize; i++)
					{
						for (int j = 0; j < TrainImgLRSize; j++)
						{
							//input1.at<float>(0, n) = (testemp.at<float>(i, j) - MinValue) / (MaxValue - MinValue);
							input1.at<float>(0, n) = testemp.at<float>(i, j);
							n++;
						}
					}
					cv::Mat testMatch = cv::Mat::ones(TrainImgHRSize + 1, (8 + 1 + 8 + 1 + 8 + 1 + 8 + 72 + 1), CV_32FC1);
					TestLum(cv::Rect(c, r, TrainImgLRSize, TrainImgLRSize)).copyTo(testemp);
					testemp.copyTo(SplitTest(cv::Rect(c*TrainImgLRSize, r*TrainImgLRSize, TrainImgLRSize, TrainImgLRSize)));
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
							testMatch.at<float>(tesr, 18 + tesc) = TestLum.at<float>(r + tesr, c + tesc);
						}
					}
					cv::Mat testemphr(TrainImgHRSize, TrainImgHRSize, CV_32FC1);
					for (int r = 0; r < TrainImgHRSize; r++)
					{
						for (int c = 0; c < TrainImgHRSize; c++)
						{
							testMatch.at<float>(0 + r, 36 + c) = TrainLum.at<float>(tempR + r, tempC + c);
							testemphr.at<float>(r, c) = TrainLum.at<float>(tempR + r, tempC + c);
						}
					}
					resize(testemphr, testemp, cv::Size(TrainImgLRSize, TrainImgLRSize), 0, 0, INTER_CUBIC);
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
					im_TrainLum[pic](cv::Range(tempR, tempR + TrainImgHRSize), cv::Range(tempC, tempC + TrainImgHRSize)).copyTo(tempHighPatch);
					float average_highpatch = Average(tempHighPatch);
					float average_lowpatch = Average(tempLowPatch);
					//resize(tempHighPatch, tempLowPatch, cv::Size(tempHighPatch.cols / DownSampleFactor, tempHighPatch.rows / DownSampleFactor), 0, 0, CV_INTER_CUBIC);
					//tempHighPatchSum = exp(param) * tempHighPatch + tempHighPatchSum;
					tempHighPatchSum = tempHighPatch - average_highpatch + average_lowpatch;
					float sum = 0.0, tempsum = 0.0;
					im_TrainLum[0](cv::Rect(tempC, tempR, TrainImgHRSize, TrainImgHRSize)).copyTo(highpatch);
					cv::resize(highpatch, lowpatch1, cv::Size(highpatch.cols / im_DownSampleFactor, highpatch.rows / im_DownSampleFactor), 0, 0, cv::INTER_AREA);
					Concatenat(lowpatch1, input1, 0);
					for (int x = 0; x<256; x++)
					{
						tempsum = abs(input1.at<float>(0, x) - TestFeatureDataOriginal.at<float>(r*(lowCols - TrainImgLRSize + 1) + c, x));
						sum = tempsum*tempsum + sum;
					}
					sum = exp(-sum / (2 * 40));
					GaussHighPatchWeight.at<float>(r, c) = sum;
				}
				//tempHighPatchSum = tempHighPatchSum / sumweight;
				/*char file_name[15];
				sprintf(file_name, "TestFinal%d.png", ttcount);
				imwrite(file_name, tempHighPatchSum*255.0);*/
				ttcount++;
				sumweight = 0, 0;
				int startIndexR = r * im_DownSampleFactor; int startIndexC = c * im_DownSampleFactor;
				for (int i = 0; i < TrainImgHRSize; i++)
				{
					for (int j = 0; j < TrainImgHRSize; j++)
					{
						//Result.at<float>(i + startIndexR, j + startIndexC) = tempHighPatchSum.at<float>(i, j)*im_HighPatchSumWeight.at<float>(i, j)*GaussHighPatchWeight.at<float>(r, c) + Result.at<float>(i + startIndexR, j + startIndexC);
						//Result.at<float>(i + startIndexR, j + startIndexC) = tempHighPatchSum.at<float>(i, j)*GaussHighPatchWeight.at<float>(r, c) + Result.at<float>(i + startIndexR, j + startIndexC);
						Result.at<float>(i + startIndexR, j + startIndexC) = tempHighPatchSum.at<float>(i, j)*im_HighPatchSumWeight.at<float>(i, j) + Result.at<float>(i + startIndexR, j + startIndexC);
						weightResult.at<float>(i + startIndexR, j + startIndexC) = im_HighPatchSumWeight.at<float>(i, j) + weightResult.at<float>(i + startIndexR, j + startIndexC);

						//weightResult.at<float>(i + startIndexR, j + startIndexC) = GaussHighPatchWeight.at<float>(r, c) + weightResult.at<float>(i + startIndexR, j + startIndexC);
					}
				}
			}
		}
		for (int i = 0; i < Result.rows; i++)
		{
			for (int j = 0; j < Result.cols; j++)
			{
				Result.at<float>(i, j) = Result.at<float>(i, j) / weightResult.at<float>(i, j);
			}
		}
		//char file_name[50];
		//sprintf(file_name, "addOriginalResult%d.png", im_PCADIMS);
		//imwrite(TEST + "addOriginalResult" + PCADIMS + ".png", Result*255.0);
		//imwrite(file_name, Result*255.0);
		im_TestLum = Result;
		cv::Mat CR = Result.clone();
		cv::Mat CB = Result.clone();
		cv::Mat CRt = cv::Mat::zeros(CR.rows, CR.cols, CV_32FC1);
		cv::Mat CBt = cv::Mat::zeros(CR.rows, CR.cols, CV_32FC1);
		resize(im_TestCR, CR, cv::Size(im_TestCB.cols * im_DownSampleFactor, im_TestCB.rows * im_DownSampleFactor), 0, 0, cv::INTER_CUBIC);
		bilateralFilter(CR, CRt, -1.0, 0.1, im_DownSampleFactor + 2);
		resize(im_TestCB, CB, cv::Size(im_TestCB.cols * im_DownSampleFactor, im_TestCB.rows * im_DownSampleFactor), 0, 0, cv::INTER_CUBIC);
		bilateralFilter(CB, CBt, -1.0, 0.1, im_DownSampleFactor + 2);
		vector< cv::Mat> LumYcYb;
		LumYcYb.push_back(Result.clone());
		LumYcYb.push_back(CRt.clone());
		LumYcYb.push_back(CBt.clone());
		cv::merge(LumYcYb, im_ResultImg);
		cv::cvtColor(im_ResultImg, im_ResultImg, cv::COLOR_YCrCb2BGR);
		im_ResultImg = im_ResultImg * 255.0;
		//imwrite(TEST + "addfinalResult.png", ResultImg);
		
		for (int m = 0; m < tcount; m++)
		{
		delete[] result[m];
		}
		delete[] result;
		LumYcYb.clear();

		//sprintf(file_name, "addfinalResult%d.png", im_PCADIMS);
		//imwrite(file_name, ResultImg);
}
void SuperOperation::ProcessImage(string ResultFile)
{
	ParameterSetting();
	ReadFilename();
	int count = 0;
	int iteration = (TESTFILENAME.size() / im_FramePersInterval)*im_FramePersInterval;
	cout << "test file numbers:"<<iteration << endl;
	//flann::Index index;
	int interal = 1;//2
	if (iteration == 1)
	{
		iteration++;
	}
	iteration = 5;
	for (int n = 0; n <= iteration; n++)
	{	
		if (n == 0)
		{
			cout << n << ":"<<endl;
			LoadTrainImg(count, interal);//2 startno number testno
			GenerateTrainSet();
			cout << "finish build training time" << endl;
			KdTreeBuild();
			cout << "finish build kdtree time" << endl;
			//PatchtableBuild();
		}
		else if (n >= count * im_FramePersInterval && n <= (count + 1) * im_FramePersInterval)
		{
			int number = n - 1;
			LoadTestImg(number);
			cout << "finish load pic" << endl;
			//cout << number << endl;
			KdTreeMethod();

			//PatchTableMethod();
			char file_name[50];
			sprintf(file_name, (ResultFile + "finalresult%d.png").c_str(), number);
			imwrite(file_name, im_ResultImg);
			cout << "write" <<number<<"pics"<< endl;
		}
		if (n == (count + 1) * im_FramePersInterval)
		{
			count++;
			//delete p;
			//p = NULL;
			/*delete im_pca;
			TrainFeatureData.release();*/
			//delete im_table;
			//im_flann_index.release();
			//im_table = NULL;
			//LoadTrainImg(count, interal);//2 startno number testno
			///cout << count << endl;

			//GenerateTrainSet();
			//KdTreeBuild();
			//PatchtableBuild();
		}
	}
	
}