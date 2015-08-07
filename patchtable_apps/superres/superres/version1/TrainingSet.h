#pragma once
#include "PatchGenerator.h"
class TrainingSet
{
public:
	TrainingSet();
	~TrainingSet();
	void Initial();
	float MatMean(Mat im);
	float Std(Mat im, float x);
	void Concatenation(Mat & KEY, Point2d *VALUE,float* HighMean, int & count);
	void Vectorize(int Row, int Col, Mat& Output,int MatSize);	float ConstrastNormal(Mat& Output);
	void FileInfo(PatchGen* LowResPatches, PatchGen* HighResPatches, float alpha, int PatchNo_row, int PatchNo_col, Mat & KEY, Point2d * VALUE,float* HighMean, int & count);
	void FileColorInfo(PatchGen* ColorLowResPatches, PatchGen* ColorHighResPatches,PatchGen* LowResPatches, PatchGen* HighResPatches, float alpha,  float colorcons,int PatchNo_row, int PatchNo_col, Mat & KEY, Point2d * VALUE, float* HighMean, int & count);
	//void AddKey_Value();//Everytime it runs, just reshape the keyMat and valueMAT.
private:
	PatchGen *m_lowResPatches, *m_highResPatches;
	PatchGen *m_colorlowResPatches, *m_colorhighResPatches;
	float m_alpha;
	int m_PatchNo_row;
	int m_PatchNo_col;
	float m_colorcons;
public:
	Mat m_Key;
	Mat *m_Value;
};

