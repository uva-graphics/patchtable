#include "stdafx.h"
#include "PatchGenerator.h"


PatchGen::PatchGen(Mat & srcImg, int size)
{
	m_srcImg = srcImg;
	//imwrite("images\\test\\patchGen.png", srcImg);
	m_patchSize = size;
}
void PatchGen::Initalize(int n)
{
	m_patch_colNo = m_srcImg.cols - m_patchSize + 1;
	m_patch_rowNo = m_srcImg.rows - m_patchSize + 1;
	m_patch = new Mat*[m_patch_rowNo];
	for (int row = 0; row < m_patch_rowNo; row++)
	{
		m_patch[row] = new Mat[m_patch_colNo];
		for (int cols = 0; cols < m_patch_colNo; cols++)
		{
			m_patch[row][cols].create(m_patchSize, m_patchSize, n);
		}
	}
	if(n==CV_32F)
	{
		for (int clipatch_r = 0; clipatch_r <m_patch_rowNo; clipatch_r++)
		{
			for (int clipatch_l = 0; clipatch_l <m_patch_colNo; clipatch_l++)
			{
				int offset_r = clipatch_r ;
				int offset_l = clipatch_l ;
				for (int patch_r = 0; patch_r < m_patchSize; patch_r++)
				{
					for (int patch_l = 0; patch_l < m_patchSize; patch_l++)
					{
						/*if (((patch_r + offset_r) == 16) && ((patch_l + offset_l) == 9))
						{
						int error = 0;
						}*/
						m_patch[clipatch_r][clipatch_l].at<float>(patch_r, patch_l) = m_srcImg.at<float>(patch_r + offset_r, patch_l + offset_l);
						/*m_patch[clipatch_r][clipatch_l].at<Vec3f>(patch_r, patch_l)[1] = m_srcImg.at<Vec3f>(patch_r + offset_r, patch_l + offset_l)[1];
						m_patch[clipatch_r][clipatch_l].at<Vec3f>(patch_r, patch_l)[2] = m_srcImg.at<Vec3f>(patch_r + offset_r, patch_l + offset_l)[2];*/
					}
				}
			}
		}
	}
	else
	{
		for (int clipatch_r = 0; clipatch_r <m_patch_rowNo; clipatch_r++)
		{
			for (int clipatch_l = 0; clipatch_l <m_patch_colNo; clipatch_l++)
			{
				int offset_r = clipatch_r ;
				int offset_l = clipatch_l ;
				for (int patch_r = 0; patch_r < m_patchSize; patch_r++)
				{
					for (int patch_l = 0; patch_l < m_patchSize; patch_l++)
					{
						m_patch[clipatch_r][clipatch_l].at<Vec3f>(patch_r, patch_l)[0] = m_srcImg.at<Vec3f>(patch_r + offset_r, patch_l + offset_l)[0];
						m_patch[clipatch_r][clipatch_l].at<Vec3f>(patch_r, patch_l)[1] = m_srcImg.at<Vec3f>(patch_r + offset_r, patch_l + offset_l)[1];
						m_patch[clipatch_r][clipatch_l].at<Vec3f>(patch_r, patch_l)[2] = m_srcImg.at<Vec3f>(patch_r + offset_r, patch_l + offset_l)[2];
					}
				}
			}
		}
	}
	//float y = m_patch[0][0].at<Vec3f>(3, 3)[0];
}

PatchGen::~PatchGen()
{
	for (int row = 0; row < m_patch_rowNo; row++)
	{
		/*for (int cols = 0; cols < m_patch_colNo; cols++)
		{
		cvReleaseMat(m_patch[row][cols]);
		}*/
		delete[]m_patch[row];
	}
	delete []m_patch;
}
Mat PatchGen::GetPatch(int row, int col)
{
	Mat result(m_patchSize, m_patchSize,CV_32F);
	result = m_patch[row][col];
	//float y = m_patch[row][col].at<Vec3f>(3, 3)[0];

	return result;
}
