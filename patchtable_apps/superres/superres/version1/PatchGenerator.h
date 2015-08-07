#pragma once
class PatchGen
{
public:
	PatchGen(Mat & srcImg, int size);
	~PatchGen();
	void Initalize(int n);
	Mat GetPatch(int row, int col);
private:
	Mat m_srcImg;
    int m_patchSize;
	
public:
	Mat **m_patch;
	int m_patch_colNo;
	int m_patch_rowNo;
};

