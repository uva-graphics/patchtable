#include "stdafx.h"
#include "TrainingSet.h"


void TrainingSet::FileInfo(PatchGen* LowResPatches, PatchGen* HighResPatches, float alpha, int PatchNo_row, int PatchNo_col, Mat & KEY, Point2d * VALUE, float* HighMean, int & count)
{
	m_lowResPatches = LowResPatches;
	m_highResPatches = HighResPatches;
	m_alpha = alpha;
	m_PatchNo_row = PatchNo_row;
	m_PatchNo_col = PatchNo_col;
	//Initial();
	Concatenation(KEY, VALUE, HighMean, count);
}
TrainingSet::TrainingSet()
{
	m_alpha = 0.0;
	m_PatchNo_row = 0;
	m_PatchNo_col = 0;
}
void TrainingSet::Vectorize(int Row, int Col, Mat& Output)
{
	Mat lowResPatch = m_lowResPatches->GetPatch(Row, Col);
	Mat highResPatch = m_highResPatches->GetPatch(Row + 1, Col + 1);
	Output = cvCreateMat(1, 58, CV_32FC1);
	int count = 0;
	for (int row = 0; row < 7; row++)
	{
		for (int col = 0; col < 7; col++)
		{
			Output.at<float>(0, count) = lowResPatch.at<float>(row, col);
			//float n0 = lowResPatch.at<Vec3f>(row, col)[channel];
			count++;
		}
	}
	for (int col = 0; col < 5; col++)
	{
		Output.at<float>(0, count) = highResPatch.at<float>(0, col) * m_alpha;
		//float lben = highResPatch.at<Vec3f>(0, col)[channel] / m_alpha;
		count++;
	}

	for (int row = 1; row < 5; row++)
	{
		Output.at<float>(0, count) = highResPatch.at<float>(row, 0) * m_alpha;
		//float lend = highResPatch.at<Vec3f>(row, 0)[channel] / m_alpha;
		count++;
	}

	//count = 0;
	//Output=Output.reshape(1);//data structure=7*7lowresolution+ first row highres+first col highres;
	/*for (int i = 0; i < 174; i++)
	{
	if (i==49||i==57||i==58||i==59)
	{
	int error = 0;
	}
	float temp = Output.at<float>(0, i);
	int q = 0;
	}*/
}
void TrainingSet::Initial()
{
	/*m_Key = Mat::zeros(m_PatchNo_row*m_PatchNo_col, 174, CV_32F);
	m_Value = new Mat[m_PatchNo_row*m_PatchNo_col];
	for (int n = 0; n < m_PatchNo_row*m_PatchNo_col; n++)
	{
	m_Value[n] = cvCreateMat(5,5,CV_32FC3);
	}*/
}

float TrainingSet::MatMean(Mat im)
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
float TrainingSet::Std(Mat im, float x)
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
	s = sqrt(s) + 0.0001;
	return s;
}
float TrainingSet::ConstrastNormal(Mat& Output)
{
	Scalar     mean;
	Scalar     stddev;
	cv::meanStdDev(Output, mean, stddev);
	float       mean_pxl = mean.val[0];
	float       stddev_pxl = stddev.val[0];
	float result = stddev_pxl + 0.0001;
	return result;
}
void TrainingSet::Concatenation(Mat & KEY, Point2d* VALUE, float* HighMean, int & count)
{

	for (int row = 0; row < m_PatchNo_row; row++)
	{
		for (int col = 0; col < m_PatchNo_col; col++)
		{

			Mat output;
			Vectorize(row, col, output);
			float x = MatMean(output);
			float mean = Std(output, x);
			//output = output / mean;
			for (int n = 0; n < 58; n++)
			{
				KEY.at<float>(count, n) = output.at<float>(0, n)/mean;
			}
			HighMean[count] = mean;
			VALUE[count] = Point2d(row, col);
			count++;

			//Mat dsttemp = KEY.row(count);             //MÎªÄ¿µÄ¾ØÕó n*m
			//output.copyTo(dsttemp);
			//float y = KEY.at<float>(0, 1);
			/*if (row == 42 && col == 193 && count == 466609)
			{
			int errpr = 0;
			Mat a = m_highResPatches->GetPatch(row + 1, col + 1) / mean;
			}*/
			/*Mat a = m_highResPatches->GetPatch(row + 1, col + 1) / mean;
			if (row==310&&col==440)
			{
			int errpr = 0;
			}
			Mat b(1, 25, CV_32FC3);
			b = VALUE(Rect(0, count, 25, 1));*/
			/*int count = 0;
			for (int row = 0; row < 5; row++)
			{
			for (int col = 0; col < 5; col++)
			{
			for (int n = 0; n < 3; n++)
			{
			b.at<Vec3f>(0, count)[n] = a.at<Vec3f>(row, col)[n];
			}
			count++;
			}
			}*/
			//VALUE.push_back(m_highResPatches->GetPatch(row + 1, col + 1) / mean);
			//VALUE(Rect(count,0,1,25))=b;
			/*if (count == 466609)
			{
			int n = 0;
			}*/
		}
	}
}

TrainingSet::~TrainingSet()
{
	//delete m_Value;
}
