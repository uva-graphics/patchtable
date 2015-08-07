#ifndef _PATCHMATCH_H_
#define _PATCHMATCH_H_

#include <opencv2/opencv.hpp>
using cv::Mat;        // Connelly: Need to use these 'using' methods instead of 'using namespace cv' to work around flann bug
using cv::Vec2f;
using cv::Vec2d;
using cv::Vec3b;
using cv::Vec3f;
using cv::Point;
using cv::Scalar;
using cv::Range;
using cv::imread;
using cv::Size;
using std::vector;
using std::string;

bool Ann(Mat Target, Mat Source, Mat& Offset, Mat& Dist, Mat Hole, Mat Guide, Mat Const, vector<vector<Vec2f>>searchSpaces, int PatchSize, int iter);
void retarget(Mat SrcImg, Mat& DstImg, Mat Cover, Mat Hole, Mat Guide, Mat Const, Mat Offset, 
	int minSize, int PatchSize, int searchWindowFacktor, int EMiterations, int ANNIterations, int coherence_spatial, int coherence_temporal, int matchIndex, float& coherenec);

#endif
