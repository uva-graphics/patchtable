#ifndef _PATCHMATCH_H_
#define _PATCHMATCH_H_

#if __APPLE__
#include <opencv2/opencv.hpp>
#else
#include <opencv.hpp>
#endif
using namespace cv;

bool Ann(Mat Target, Mat Source, Mat& Offset, Mat& Dist, Mat Hole, Mat Guide, Mat Const, vector<vector<Vec2f> >searchSpaces, int PatchSize, int iter);
void retarget(Mat SrcImg, Mat& DstImg, Mat Hole, Mat Guide, Mat Const, Mat Offset, int minSize, int PatchSize, int searchWindowFacktor, int EMiterations, int ANNIterations);

#endif
