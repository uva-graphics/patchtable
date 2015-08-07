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
	void Vectorize(int Row, int Col, Mat& Output);	float ConstrastNormal(Mat& Output);
	void FileInfo(PatchGen* LowResPatches, PatchGen* HighResPatches, float alpha, int PatchNo_row, int PatchNo_col, Mat & KEY, Point2d * VALUE,float* HighMean, int & count);
	//void AddKey_Value();//Everytime it runs, just reshape the keyMat and valueMAT.
private:
	PatchGen *m_lowResPatches, *m_highResPatches;
	float m_alpha;
	int m_PatchNo_row;
	int m_PatchNo_col;
public:
	Mat m_Key;
	Mat *m_Value;
};

