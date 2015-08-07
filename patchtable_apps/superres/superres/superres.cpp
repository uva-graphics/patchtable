// SupeRes.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "PatchGenerator.h"
#include <iostream>
#include <fstream>
#include "vld.h"
using namespace std;
#include "TrainingSet.h"
vector<string> TRAINFILENAME;
//string INPUTFILNAME = "images\\input\\input.png";
//string INPUTFILNAME = "images\\input\\123.png";
string INPUTFILNAME = "images\\input\\100075.jpg";
string TEST = "images\\test\\";
string KEY = "E:\\KEY.txt";
string OUTPUT = "E:\\OUTPUT.txt";
float ALPHA = 0.54;
Mat Key;
//vector<Mat> Values;
//Mat Values;
Mat Result;
Mat *InputImg, *InputHigh;
Point2d* Valuse;
int* fileSum; float * HighMean;
int traiNo;
void WriteFile(int rows, int cols, float **data, string NAME){
	/*int m = 10;
	int n = 20;
	*/
	//float **data = new float*[rows];
	/*for (int i = 0; i <rows; i++){
	data[i] = new float[cols];
	for (int j = 0; j < cols; j++){
	data[i][j] = 1.1f;
	}
	}
	*/
	float x_min = 500.0, x_max = -500.0;
	ofstream out(NAME);
	for (int i = 0; i < rows; i++){
		for (int j = 0; j < cols - 1; j++){
			out << data[i][j] << " ";
			if (data[i][j]>x_max)
			{
				x_max = data[i][j];
			}
			if (data[i][j]<x_min)
			{
				x_min = data[i][j];
			}
		}
		out << data[i][cols - 1] << endl;
		if (data[i][cols - 1]>x_max)
		{
			x_max = data[i][cols - 1];
		}
		if (data[i][cols - 1] < x_min)
		{
			x_min = data[i][cols - 1];
		}
	}
	out << x_min << " ====" << x_max;
	out.close();
}
int LoadTrainPic()
{
	/*TRAINFILENAME.push_back("images\\train\\flowers1.png");
	TRAINFILENAME.push_back("images\\train\\flowers2.png");
	TRAINFILENAME.push_back("images\\train\\flowers3.png");*/
	TRAINFILENAME.push_back("images\\train\\15004.jpg");
	TRAINFILENAME.push_back("images\\train\\22093.jpg");
	TRAINFILENAME.push_back("images\\train\\100080.jpg");
	TRAINFILENAME.push_back("images\\train\\105053.jpg");
	TRAINFILENAME.push_back("images\\train\\108041.jpg");
	TRAINFILENAME.push_back("images\\train\\113009.jpg");
	TRAINFILENAME.push_back("images\\train\\113016.jpg");
	TRAINFILENAME.push_back("images\\train\\159091.jpg");
	TRAINFILENAME.push_back("images\\train\\227040.jpg");
	TRAINFILENAME.push_back("images\\train\\247085.jpg");
	/*TRAINFILENAME.push_back("images\\train\\274007.jpg");
	TRAINFILENAME.push_back("images\\train\\277095.jpg");
	TRAINFILENAME.push_back("images\\train\\302003.jpg");
	TRAINFILENAME.push_back("images\\train\\317080.jpg");
	TRAINFILENAME.push_back("images\\train\\365025.jpg");
	TRAINFILENAME.push_back("images\\train\\368016.jpg");*/

	traiNo = TRAINFILENAME.size();
	InputImg = new Mat[traiNo];
	InputHigh = new Mat[traiNo];
	/*TRAINFILENAME.push_back("images\\train\\23025.jpg");
	TRAINFILENAME.push_back("images\\train\\24004.jpg");
	TRAINFILENAME.push_back("images\\train\\25098.jpg");
	TRAINFILENAME.push_back("images\\train\\27059.jpg");*/
	int num = 0;
	fileSum = new int[traiNo + 1];
	fileSum[0] = -1;
	for (int i = 0; i < traiNo; i++)
	{
		InputImg[i] = imread(TRAINFILENAME[i], CV_LOAD_IMAGE_GRAYSCALE);
		int m_patch_colNo = InputImg[i].cols - 7 + 1;
		int m_patch_rowNo = InputImg[i].rows - 7 + 1;
		num = m_patch_colNo*m_patch_rowNo + num;
		fileSum[i + 1] = num;
		//InputHigh[i].create(InputImg[i].rows,InputImg[i].cols,);
	}
	//delete InputImg;
	return num;

}
void GeneraLHImg(Mat & LowResImg, Mat & HighResImg, Mat &SrcImg, Mat &cubicImg)
{
	int cols = SrcImg.cols;
	int rows = SrcImg.rows;
	Mat blurImg, cubicImgTemp, Difference, blur, blurlow, blurtemp, blurlowlow, lowlow;
	GaussianBlur(SrcImg, blur, Size(7, 7), 1, 0.0, BORDER_REPLICATE);
	//imwrite(TEST + "00blur.png", blur);
	//imwrite(TEST + "0blur.png", blur);
	resize(blur, blurImg, Size(blur.cols / 4, blur.rows / 4), 0, 0, CV_INTER_CUBIC);
	/*float r = blurImg.at<Vec3f>(50, 45)[0];
	float g = blurImg.at<Vec3f>(50, 45)[1];
	float b = blurImg.at<Vec3f>(50, 45)[2];*/
	resize(blurImg, cubicImg, Size(blurImg.cols * 4, blurImg.rows * 4), 0, 0, CV_INTER_CUBIC);//WHY /4 CAN'T WORK??

	//imwrite(TEST + "traincubic.png", cubicImg);
	HighResImg = SrcImg - cubicImg;
	//imwrite(TEST + "highpass.png", HighResImg);

	
	//pyrDown(SrcImg, blurlow, Size(SrcImg.cols / 2, SrcImg.rows / 2));
	//resize(blurlow, cubicImgTemp, Size(blurlow.cols * 2, blurlow.rows * 2), 0, 0, CV_INTER_CUBIC);//WHY /4 CAN'T WORK??
	GaussianBlur(SrcImg, blurtemp, Size(25, 25), 5, 0.0, BORDER_REPLICATE);
	//pyrDown(blurtemp, blurlowlow, Size(blurtemp.cols / 2, blurtemp.rows / 2));
	//	imwrite(TEST + "1sizedown.png", blurImg);
	//resize(blurlowlow, lowlow, Size(blurlowlow.cols * 2, blurlowlow.rows * 2), 0, 0, CV_INTER_CUBIC);//WHY /4 CAN'T WORK??
	//imwrite(TEST + "4blur.png", cubicImgTemp);
	LowResImg = cubicImg - blurtemp;
	cubicImg = blurtemp;
	//imwrite(TEST + "bandpass.png", LowResImg);
	/*imwrite(TEST + "5L+blur.png", LowResImg + cubicImgTemp);
	imwrite(TEST + "6H+cubic.png", HighResImg + cubicImg);*/
}
void GeneraTrainSet()
{
	int num = LoadTrainPic();
	int traiNo = TRAINFILENAME.size();
	TrainingSet *trainSet = new TrainingSet();
	int keycount = 0, currentImgNo = 0;
	//Values(num);
	Key = cvCreateMat(num, 58, CV_32F);
	HighMean = new float[num];
	Valuse = new Point2d[num];
	//Values = cvCreateMat(num, 25, CV_32FC3);
	for (int i = 0; i < traiNo; i++)
	{
		//Mat srcImg = imread(TRAINFILENAME[i], 1);
		//Mat srcImg = imread("images\\input\\123.jpg",1);//CV_32S
		//Mat srcrgb = imread("images\\input\\123.jpg", CV_LOAD_IMAGE_GRAYSCALE);//CV_8U
		Mat srcImg = InputImg[i];
		srcImg.convertTo(srcImg, CV_32FC1);
		for (int row = 0; row < srcImg.rows; row++)
		{
			for (int col = 0; col < srcImg.cols; col++)
			{
				srcImg.at<float>(row, col) = srcImg.at<float>(row, col) / 255.0;
			}
		}
		Mat lowResImg, highResImg, cubicImg;
		GeneraLHImg(lowResImg, highResImg, srcImg, cubicImg);
		char temps[20];
		sprintf(temps, "%d_lowresImg.png", i);
		string s(temps);
		imwrite(TEST + s, lowResImg*255.0);
		sprintf(temps, "%d_highImg.png", i);
		string c(temps);
		imwrite(TEST + c, highResImg*255.0);

		PatchGen *lowResPatches = new PatchGen(lowResImg, 7);
		lowResPatches->Initalize();
		PatchGen *higResPatches = new PatchGen(highResImg, 5);
		higResPatches->Initalize();
		InputHigh[i] = highResImg.clone();
		/*if (i == 4){
			int error = 1;
			}*/
		trainSet->FileInfo(lowResPatches, higResPatches, ALPHA, lowResPatches->m_patch_rowNo, lowResPatches->m_patch_colNo, Key, Valuse, HighMean, keycount);
		delete lowResPatches;
		delete higResPatches;
	}
	delete trainSet;
}
int InWhichPic(int index)
{
	int n = 0;
	for (int x = 0; x < traiNo; x++)
	{
		if (index >= fileSum[x] && index < fileSum[x + 1])
		{
			return x;
		}

	}
	return -1;

}
float MatMean(Mat im)
{
	int n = 0;
	float s = 0.0;
	for (int r = 0; r < im.rows; r++)
	{
		for (int c = 0; c < im.cols; c++)
		{
			s = im.at<float>(r, c) + s;
			n++;
		}
	}
	s = s / n;
	return s;
}
float Std(Mat im, float x)
{
	float s = 0.0; int n = 0;
	for (int r = 0; r < im.rows; r++)
	{
		for (int c = 0; c < im.cols; c++)
		{
			s = (im.at<float>(r, c) - x)*(im.at<float>(r, c) - x);
			n++;
		}
	}
	s = s / n;
	s = sqrt(s)+0.0001;
	return s;
}
Mat  IMconstrast(Mat &im, int wsize)
{
	int height = im.rows, width = im.cols;

	Mat H;
	H = Mat::zeros(height, width, CV_32F);
	for (int r = 0; r < height; r++)
	{
		for (int c = 0; c < width; c++)
		{
			int x1 = max(c - wsize, 0);
			int y1 = max(r - wsize, 0);
			int x2 = min(c + wsize, width - 1);
			int y2 = min(r + wsize, height - 1);
			Mat patch = im(Rect(x1, y1, x2 - x1, y2 - y1));
			float mean = MatMean(patch);
			float x = Std(patch, mean);
			H.at<float>(r, c) = Std(patch, mean);
		}
	}
	/*double minVal;
	double maxVal;
	Point minLoc;
	Point maxLoc;

	minMaxLoc(H, &minVal, &maxVal, &minLoc, &maxLoc);

	cout << "min val : " << minVal << endl;
	cout << "max val: " << maxVal << endl;*/
	for (int r = 0; r < height; r++)
	{
		for (int c = 0; c < width; c++)
		{
			H.at<float>(r, c) = im.at<float>(r, c) / (H.at<float>(r, c) + 0.01);
		}
	}
	return H;
}
void SuperRe()
{
	//Key = cvCreateMat(446508, 174, CV_32F);
	/*******************************************/
	TickMeter tm;
	tm.start();
	float** data;
	vector<Mat> CrCb;
	//data = new float*[Key.rows];
	/*for (int i = 0; i < Key.rows; i++){
		data[i] = new float[Key.cols];
		for (int j = 0; j < Key.cols; j++){
		data[i][j] = -500.0;
		}
		}
		for (int row = 0; row < Key.rows; row++)
		{
		for (int col = 0; col < Key.cols; col++)
		{
		data[row][col] = Key.at<float>(row, col);
		}
		}
		WriteFile(Key.rows, Key.cols, data, KEY);
		for (int i = 0; i < Key.rows; i++){
		delete[]data[i];
		}
		delete[]data; */
	/*******************************************/
	//ifstream in(KEY);
	//char line[10240];
	//in.getline(line, 10240);
	//cout << line;
	//int lineno = 0;
	//while (in.getline(line, 10240)){
	//	stringstream ss2;
	//	ss2 << line;
	//	for (int j = 0; j < 174; j++){
	//		ss2 >> Key.at<float>(lineno, j);// data[lineno][j];
	//	}
	//	lineno++;
	//}
	/*for (int row = 0; row < Key.rows; row++)
	{
	for (int col = 0; col < Key.cols; col++)
	{
	Key.at<float>(row, col) = data[row][col];
	}
	}*/
	Mat inputImg = imread(INPUTFILNAME, 1);
	inputImg.convertTo(inputImg, CV_32FC3);
	for (int row = 0; row < inputImg.rows; row++)
	{
		for (int col = 0; col < inputImg.cols; col++)
		{
			inputImg.at<Vec3f>(row, col)[0] = inputImg.at<Vec3f>(row, col)[0]/255.0;
			inputImg.at<Vec3f>(row, col)[1] = inputImg.at<Vec3f>(row, col)[1] / 255.0;
			inputImg.at<Vec3f>(row, col)[2] = inputImg.at<Vec3f>(row, col)[2] / 255.0;

		}
	}
	Mat inputImage, Input, rycImg, temp, Cr, Cb, srcOfMerge, finalMerge;
	GaussianBlur(inputImg, inputImage, Size(7, 7), 1, 0.0, BORDER_REPLICATE);
	resize(inputImage, Input, Size(inputImage.cols / 4, inputImage.rows / 4), 0, 0, CV_INTER_CUBIC);
	imwrite(TEST + "inputImage.png", Input*255.0);
	
	cvtColor(inputImg, rycImg, CV_BGR2YCrCb);
	split(rycImg, CrCb);
	Mat SrcImg = CrCb[0];
	GaussianBlur(CrCb[1], CrCb[1], Size(7, 7), 1, 0.0, BORDER_REPLICATE);
	resize(CrCb[1], temp, Size(inputImage.cols / 4, inputImage.rows / 4), 0, 0, CV_INTER_CUBIC);
	resize(temp, CrCb[1], Size(inputImage.cols, inputImage.rows), 0, 0, CV_INTER_CUBIC);

	GaussianBlur(CrCb[2], CrCb[2], Size(7, 7), 1, 0.0, BORDER_REPLICATE);
	resize(CrCb[2], temp, Size(inputImage.cols / 4, inputImage.rows / 4), 0, 0, CV_INTER_CUBIC);
	resize(temp, CrCb[2], Size(inputImage.cols, inputImage.rows), 0, 0, CV_INTER_CUBIC);

	Mat lowResImg, highResImg, cubicImg, bandpass, AddDifference, blurImg;
	GeneraLHImg(bandpass, highResImg, SrcImg, cubicImg);
	lowResImg = bandpass + cubicImg;
	imwrite(TEST + "inputbandpass.png", bandpass);
	imwrite(TEST + "inputbandpass1.png", bandpass*255.0);

	CrCb[0] = lowResImg;
	merge(CrCb, srcOfMerge);
	cvtColor(srcOfMerge, finalMerge, CV_YCrCb2BGR);
	imwrite(TEST + "inputlowpasscompa.png", finalMerge*255.0);
	Mat foo = IMconstrast(bandpass, 3);
	imwrite(TEST + "inputbandpasscn1.png", foo*255.0);
	Mat abss = abs(foo);
	Mat contrast(abss.rows, abss.cols, CV_32F);
	float max = 0;
	for (int r = 0; r < abss.rows; r++)
	{
		for (int c = 0; c < abss.cols; c++)
		{
			if (abss.at<float>(r, c)>max)
			{
				max = abss.at<float>(r, c);
			}
		}
	}

	for (int r = 0; r < abss.rows; r++)
	{
		for (int c = 0; c < abss.cols; c++)
		{
			float x = foo.at<float>(r, c);
			float y = foo.at<float>(r, c) / (max + 0.5);

			contrast.at<float>(r, c) = foo.at<float>(r, c) / (max + 0.5) + 0.5;
			float z = contrast.at<float>(r, c);
			float xx = 0;
		}
	}
	imwrite(TEST + "inputbandpassContrast.png", contrast*225.0);
	//Result.create(bandpass.rows, bandpass.cols, CV_32F);// Scalar::all(128));// ::zeros(lowResImg.rows, lowResImg.cols, CV_32FC3);R
	//Result = bandpass.clone();// setTo(Scalar::all(0));
	//Result = highResImg.clone();
	Result = highResImg.clone();
	Result.setTo(Scalar::all(0));
	//Mat subResult = Result.clone();
	cv::flann::Index kdtree(Key, cv::flann::KDTreeIndexParams(1));

	int iterate_r = Result.rows - 6, iterate_c = Result.cols - 6;
	vector<int> indices(1);
	vector<float> dists(1);
	vector<int>::iterator it;
	map<int, float>::iterator it_map;
	map<int, float> vectorList;
	ofstream fout("test.txt");
	if (!fout)
	{
		cout << "File Not Opened" << endl;  return;
	}
	int index = 0;

	/**************************************/
	/*ifstream in("E:\\testyy.txt");
	int** data;
	char line[10240];
	in.getline(line, 10240);*/
	/*int datanum = 0;
	for (int r = 0; r < iterate_r; r = r + 4)
	{
	for (int c = 0; c < iterate_c; c = c + 4)
	{
	datanum++;
	}
	}*/
	/**************************************/
	//int rowr = 0;
	////float** data;
	//data = new float*[datanum];
	//for (int i = 0; i < datanum; i++){
	//	data[i] = new float[174];
	//	for (int j = 0; j < 174; j++){
	//		data[i][j] = -500.0;
	//	}
	//}
	/**************************************/
	for (int r = 0; r < iterate_r; r = r + 4)
	{
		/******************************/
		//in.getline(line, 10240);
		//stringstream ss2;
		//ss2 << line;
		/******************************/
		for (int c = 0; c < iterate_c; c = c + 4)
		{
			Mat lowres = bandpass(Rect(c, r, 7, 7));
			Mat higres = Result(Rect(c + 1, r + 1, 5, 5));
			//Mat subhigres = subResult(Rect(c + 1, r + 1, 5, 5));
			Mat Output = cvCreateMat(1, 58, CV_32FC1);
			int count = 0;
			/*for (int channel = 0; channel < 3; channel++)
			{*/
			for (int row = 0; row < 7; row++)
			{
				for (int col = 0; col < 7; col++)
				{
					Output.at<float>(0, count) = lowres.at<float>(row, col);
					count++;
				}
			}
			for (int col = 0; col < 5; col++)
			{
				Output.at<float>(0, count) = higres.at<float>(0, col) * ALPHA;
				count++;
			}

			for (int row = 1; row < 5; row++)
			{
				Output.at<float>(0, count) = higres.at<float>(row, 0) * ALPHA;
				count++;
			}
			/*}*/
			float x = MatMean(Output);
			float result = Std(Output, x);
			Output = Output / result;
			/*for (int n = 0; n < 174; n++)
			{
			data[rowr][n] = Output.at<float>(0, n);
			}
			rowr++;*/
			kdtree.knnSearch(Output, indices, dists, 1, cv::flann::SearchParams(32));

			it = indices.begin();
			index = *it;
			int pic = InWhichPic(index);
			if (pic == -1 || pic > traiNo)
			{
				int error = 0;
			}
			Point2d point = Valuse[index];
			int rh = -1, ch = 0;
			for (int patch_r = point.x + 1; patch_r < point.x + 6; patch_r++)
			{
				rh++; ch = 0;
				for (int patch_l = point.y + 1; patch_l < point.y + 6; patch_l++)
				{
					higres.at<float>(rh, ch) = ((InputHigh[pic].at<float>(patch_r, patch_l)) / HighMean[index])* result;
					/*higres.at<Vec3f>(rh, ch)[1] = ((InputHigh[pic].at<Vec3f>(patch_r, patch_l)[1]) / HighMean[index])* result;
					higres.at<Vec3f>(rh, ch)[2] = ((InputHigh[pic].at<Vec3f>(patch_r, patch_l)[2]) / HighMean[index])* result;*/
					ch++;
				}
			}
			//Mat temp=
			/******************/
			//ss2 >> index;
			/******************/
			vectorList.insert(make_pair(index, result));
			//higres = Values[index] * result;
			//Mat temp = Values(Rect(index, 0, 1, 25)); int cc = 0;
			/*for (int row = 0; row < 5; row++)
			{
			for (int col = 0; col < 5; col++)
			{
			for (int n = 0; n < 3; n++)
			{
			higres.at<Vec3f>(row, col)[n] = temp.at<Vec3f>(0, cc)[n] * result;
			}
			cc++;
			}
			}*/
			//higres = Values(Rect(index, 0, 1, 25))*result;
			//if (r<20&&c<21)
			//{
			char temps[20];
			sprintf(temps, "%d ", index);
			string s(temps);
			fout << s << " ";
			/*}*/

		}
		/*if (r < 20)
		{*/
		fout << endl;
		/*}*/
	}
	/*************************************/
	/*WriteFile(3025, 174, data,OUTPUT);
	for (int i = 0; i <3025; i++){
	delete[]data[i];
	}
	delete[]data;*/
	/*************************************/
	int part1 = 0, part2 = 0, part3 = 0, part4 = 0;
	int *sumtemp = new int[traiNo];
	for (int n = 0; n < traiNo; n++)
	{
		sumtemp[n] = 0;
	}
	for (it_map = vectorList.begin(); it_map != vectorList.end(); it_map++)
	{
		int index = it_map->first;
		int pic = InWhichPic(index);
		//	//if (index > 0 && index < 148836)
		sumtemp[pic]++;
		//	char temps[64];
		//	sprintf(temps, "[%d]%d.png", pic, it_map->first);
		//	string s(temps);
		//	Mat higres(5, 5, CV_32FC3);
		//	Point2d point = Valuse[index]; int rr = -1, cr = 0;
		//	for (int patch_r = point.x + 1; patch_r <point.x+ 6; patch_r++)
		//	{
		//		rr++; cr = 0;
		//		for (int patch_l = point.y + 1; patch_l < point.y+6; patch_l++)
		//		{
		//			higres.at<Vec3f>(rr, cr)[0] = InputHigh[pic].at<Vec3f>(patch_r, patch_l)[0];
		//			higres.at<Vec3f>(rr, cr)[1] = InputHigh[pic].at<Vec3f>(patch_r, patch_l)[1];
		//			higres.at<Vec3f>(rr, cr)[2] = InputHigh[pic].at<Vec3f>(patch_r, patch_l)[2];
		//			float n = HighMean[index];
		//			cr++;
		//		}
		//	}
		//	imwrite("images\\test\\input" + s, higres);
	}
	for (int n = 0; n < traiNo; n++)
	{
		fout << sumtemp[n] << endl;
	}
	fout.close();
	Mat final, finalresult;
	imwrite(TEST + "addOriginalResult.png", Result*255.0);
	add(lowResImg, Result, Result);
	imwrite(TEST + "addResult.png", Result*255.0);
	CrCb[0] = Result;
	merge(CrCb, finalresult);
	cvtColor(finalresult, finalMerge, CV_YCrCb2BGR);
	imwrite(TEST + "addfinalResult.png", finalMerge*255.0);
	tm.stop();
	cout << "count=" << tm.getCounter() << ",process time=" << tm.getTimeSec() << endl;
	Sleep(5000);
	CrCb.clear();
	vectorList.clear();
	/*cvMerge(Y, Cr, Cb, NULL, frame)
		cvtColor(frame, frame, CV_YCrCb2BGR)*/
}
int _tmain(int argc, _TCHAR* argv[])
{
	GeneraTrainSet();
	SuperRe();
	return 0;
}
