#include "reshuffle.h"
#include "../../patchtable/patchtable.h"
#include<ctime>
#define MAX_DIST 1000000 //maximum distance for the distance function
#define ZERO 0.00000000001
//#define _DEBUG_ //show all src, dst, hole, offset, guide, constraint maps.
//#define _COMPARE_UP_LAYER_ //compare the patchtable result with the uplayer offset map
//#define _FIX_ //fix the coherence
#define _PATCHTABLE_//use patchtable, or use patchmatch
//for debug
void showGUI(Mat ori, Mat dst){
	int row = ori.rows, col = ori.cols;
	Mat dstImg = Mat(row, col, CV_8UC3);
	for (int i = 0; i < row; i++){
		for (int j = 0; j < col; j++){
			dstImg.at<Vec3b>(i, j) = dst.at<Vec3f>(i*dst.rows / row, j*dst.cols / col)*255.f;
		}
	}
	imshow("GUI", dstImg);
	cvWaitKey(1);
}

void showAll(Mat ori, Mat src, Mat dst, Mat hole, Mat guide, Mat off, Mat con){
#ifdef _DEBUG_
	int row = ori.rows, col = ori.cols;
	Mat srcImg = Mat(row, col, CV_8UC3);
	Mat dstImg = Mat(row, col, CV_8UC3);
	Mat holeImg = Mat(row, col, CV_8UC3);
	Mat guideImg = Mat(row, col, CV_8UC3);
	Mat offImg = Mat(row, col, CV_8UC3);
	Mat conImg = Mat(row, col, CV_8UC3);
	for (int i = 0; i < row; i++){
		for (int j = 0; j < col; j++){
			srcImg.at<Vec3b>(i, j) = src.at<Vec3f>(i*src.rows / row, j*src.cols / col)*255.f;
			dstImg.at<Vec3b>(i, j) = dst.at<Vec3f>(i*dst.rows / row, j*dst.cols / col)*255.f;
			if (hole.at<uchar>(i*hole.rows / row, j*hole.cols / col)){
				holeImg.at<Vec3b>(i, j) = Vec3b(0, 0, 0);
			}
			else{
				holeImg.at<Vec3b>(i, j) = Vec3b(255, 255, 255);
			}
			if (guide.at<uchar>(i*guide.rows / row, j*guide.cols / col)){
				guideImg.at<Vec3b>(i, j) = Vec3b(0, 0, 0);
			}
			else{
				guideImg.at<Vec3b>(i, j) = Vec3b(255, 255, 255);
			}
			if (con.at<uchar>(i*con.rows / row, j*con.cols / col)){
				conImg.at<Vec3b>(i, j) = Vec3b(0, 0, 0);
			}
			else{
				conImg.at<Vec3b>(i, j) = Vec3b(255, 255, 255);
			}
			int r = (off.at<Vec2f>(i*off.rows / row, j*off.cols / col)[0] + (src.cols + dst.cols)) * 255 / (2 * (src.cols + dst.cols));
			int g = (off.at<Vec2f>(i*off.rows / row, j*off.cols / col)[1] + (src.rows + dst.rows)) * 255 / (2 * (src.rows + dst.rows));
			offImg.at<Vec3b>(i, j) = Vec3b(0, g, r);
		}
	}
	imshow("src", srcImg);
	imshow("dst", dstImg);
	imshow("hole", holeImg);
	imshow("guide", guideImg);
	imshow("off", offImg);
	imshow("const", conImg);
	cvWaitKey();
#endif
}


//if the patches not exsit or not equal size, this function will return -1
double distance(Mat Target, Mat Source, int tx, int ty, int sx, int sy, int PatchSize, double abort = MAX_DIST);
void pyrDownHole(Mat HoleIn, Mat GuideIn, Mat& HoleOut, Mat& GuideOut);
void pyrDownOffs(Mat OffIn, Mat ConstIn, Mat& OffOut, Mat& ConstOut);
void pyrUpOffset(Mat OffsetIn, Mat& OffsetOut, Mat HoleOut, Mat GuideOut, Mat ConstOut, vector<vector<Vec2f>> searchSpaces);
void generateImage(Mat Offset, Mat SrcImg, Mat& DstImg, int patchSize, bool up = false);

PatchTable<float, float>* build_table(Mat B, Mat Hole){
	PatchTableParams *p = new PatchTableParams();
	p->calc_exact_dist = true;
	p->set_speed(9);
	p->coherence_spatial=4.0;
	Mat B_;
	Array<int32_t> allow(Hole.rows, Hole.cols);
	for (int row = 0; row < Hole.rows; row++){
		for (int col = 0; col < Hole.cols; col++){
			if (!Hole.at<uchar>(row, col)){
				allow(row, col) = 1;
			}
			else{
				allow(row, col) = 0;
			}
		}
	}
	B.convertTo(B_, CV_8UC3, 255.f);
	imwrite("B.png", B_);
	Array<float> b(load_color_image<float>("B.png"));
	PatchTable<float, float>* table = new PatchTable<float, float>(p, b, &allow);
	return table;
}

void table_match(PatchTable<float, float> *table, Mat A, Mat B, Mat& offset, Mat Hole, int PatchSize){
	Array<double> ann;
	Mat A_;
	A.convertTo(A_, CV_8UC3, 255.f);
	imwrite("A.png", A_);
	Array<float> a(load_color_image<float>("A.png"));
	for(int i=0;i<1;i++){
		table->lookup(a, ann, 0, &ann);
	}
	
	for (int i = 0; i < offset.rows - 7; i++){
		for (int j = 0; j < offset.cols - 7; j++){
			int tx = j + 3, ty = i + 3;
			Vec2f off = offset.at<Vec2f>(ty, tx);
			float x = float(ann(i, j, NNF_X));
			float y = float(ann(i, j, NNF_Y));
#ifdef _COMPARE_UP_LAYER_
			int sx = tx + off[0], sy = ty + off[1];
			double ori_dist = distance(A, B, tx, ty, sx, sy, PatchSize);
			double dist = distance(A, B, tx, ty, x + 3, y + 3, PatchSize);
			if (!Hole.at<uchar>(sy, sx) && dist < ori_dist*0.5){
				offset.at<Vec2f>(ty, tx) = Vec2f(x - j, y - i);
			}
#else
			offset.at<Vec2f>(ty, tx) = Vec2f(x - j, y - i);
#endif
		}
	}
#ifdef _FIX_
	for (int i = 4; i < offset.rows - 5; i++){
		for (int j = 4; j < offset.cols - 5; j++){
			Vec2f around[8];
			around[0] = offset.at<Vec2f>(i - 1, j - 1);
			around[1] = offset.at<Vec2f>(i, j - 1);
			around[2] = offset.at<Vec2f>(i + 1, j - 1);
			around[3] = offset.at<Vec2f>(i - 1, j);
			around[4] = offset.at<Vec2f>(i + 1, j);
			around[5] = offset.at<Vec2f>(i - 1, j + 1);
			around[6] = offset.at<Vec2f>(i, j + 1);
			around[7] = offset.at<Vec2f>(i + 1, j + 1);
			int index = 0;
			int num[8] = { 0 };
			for (int m = 0; m < 8; m++){
				bool found = false;
				for (int n = 0; n < index; n++){
					if (around[m] == around[n]){
						num[n]++;
						found = true;
						break;
					}
				}
				if (!found){
					around[index] = around[m];
					num[index] = 1;
					index++;
				}
			}
			int max = 0;
			bool valid = false;
			Vec2f centre = offset.at<Vec2f>(i, j);
			for (int k = 0; k < index; k++){
				if (around[k] == centre && num[k] > 1){
					valid = true;
					break;
				}
				if (num[k] > max){
					max = k;
				}
			}
			if (!valid){
				offset.at<Vec2f>(i, j) = around[max];
			}
		}
	}
#endif
}
//void table_match(Mat A, Mat B, Mat& offset){
//	PatchTableParams *p = new PatchTableParams();
//	p->calc_exact_dist = true;
//	Mat A_,B_;
//	A.convertTo(A_, CV_8UC3, 255.f);
//	B.convertTo(B_, CV_8UC3, 255.f);
//	imwrite("A.png", A_);
//	imwrite("B.png", B_);
//	Array<float> a(load_color_image<float>("A.png"));
//	Array<float> b(load_color_image<float>("B.png"));
//	PatchTable<float, float> table(p, b);
//	Array<double> ann;
//	double latest_time = 1e100;
//	for (int iter = 0; iter < 2; iter++){
//		latest_time = table.lookup(a, ann);
//	}
//	save_color_image(ann, "out.png");
//	for (int i = 0; i < offset.rows-7; i++){
//		for (int j = 0; j < offset.cols-7; j++){
//			float x = float(ann(i, j, NNF_X));
//			float y = float(ann(i, j, NNF_Y));
//			//printf("x:%d,y:%d\n", x, y);
//			offset.at<Vec2f>(i + 3, j + 3)[0] = x - j;
//			offset.at<Vec2f>(i + 3, j + 3)[1] = y - i;
//		}
//	}
//}

bool Ann(Mat Target, Mat Source, Mat& Offset, Mat& Dist, Mat Hole, Mat Guide, Mat Const, vector<vector<Vec2f>>searchSpaces, int PatchSize, float searchWindow, float alpha, int iter){
	//check format
	if (Target.type() != CV_32FC3){
		fprintf(stderr, "Target type not support.%d,%d\n", Target.type(), CV_32FC3);
		return false;
	}
	if (Source.type() != CV_32FC3){
		fprintf(stderr, "Source type not support.\n");
		return false;
	}
	if (Offset.type() != CV_32FC2){
		fprintf(stderr, "Offset type not support.\n");
		return false;
	}
	//decide loop order
	int xfirst, xlast, yfirst, ylast, step;
	if (iter % 2 == 0){
		xfirst = 0, xlast = Offset.cols;
		yfirst = 0, ylast = Offset.rows;
		step = 1;
	}
	else{
		xfirst = Offset.cols - 1, xlast = -1;
		yfirst = Offset.rows - 1, ylast = -1;
		step = -1;
	}
	//begin the loop
	for (int tx = xfirst; tx != xlast; tx += step){
		for (int ty = yfirst; ty != ylast; ty += step){
			//pass constraint patch
			if (Const.at<uchar>(ty, tx)){
				continue;
			}
			//get current distance
			int guide = Guide.at<uchar>(ty, tx);
			if (guide >= searchSpaces.size()){
				guide = 0;
			}
			Vec2f offset = Offset.at<Vec2f>(ty, tx);
			int sx = offset[0] + tx, sy = offset[1] + ty;
			double bestDistance = Dist.at<double>(ty, tx);
			if (bestDistance < 0){
				bestDistance = distance(Target, Source, tx, ty, sx, sy, PatchSize);
			}
			if (bestDistance < 0){
				//				fprintf(stderr,"Patch T(%d,%d)->S(%d,%d) are not same size.\n",tx,ty,sx,sy);
				if (searchSpaces[guide].size() == 0){
					fprintf(stderr, "Patch T(%d,%d) belongs to Guide[%d] (Guide[%d] is empty).\n", tx, ty, guide, guide);
					guide = 0;
				}
				if (searchSpaces[0].size() == 0){
					fprintf(stderr, "Whole Target is hole,can't match anything.\n");
					return false;
				}
				//choose a valid init for this Target patch
				bool valid = false;
				for (int i = searchSpaces[guide].size(); i > 0; i--){
					int index = rand() % i;
					Vec2f s = searchSpaces[guide][index];
					sx = s[0], sy = s[1];
					bestDistance = distance(Target, Source, tx, ty, sx, sy, PatchSize);
					if (bestDistance >= 0){
						Offset.at<Vec2f>(ty, tx) = s - Vec2f(tx, ty);
						valid = true;
						break;
					}
					Vec2f temp = searchSpaces[guide][index];
					searchSpaces[guide][index] = searchSpaces[guide][i - 1];
					searchSpaces[guide][i - 1] = temp;
				}
				if (valid == false){
					fprintf(stderr, "Patch T(%d,%d) have no valid match.\n", tx, ty);
				}
			}
			//convergence
			if (abs(bestDistance) < ZERO){
				Dist.at<double>(ty, tx) = bestDistance;
				continue;
			}

			//propergate
			if (tx != xfirst){
				offset = Offset.at<Vec2f>(ty, tx - step);
				sx = tx + offset[0], sy = ty + offset[1];
				if (sx<0 || sx>Hole.cols - 1){
					sx = sx - step;
					if (sx<0 || sx>Hole.cols - 1 || sy<0 || sy>Hole.rows - 1){
						fprintf(stderr, "!!!!!!!!!!!!!!%d,%d->%d,%d\n", tx, ty, sx, sy);
					}

					if (!Hole.at<uchar>(sy, sx) && (guide == 0 || Guide.at<uchar>(sy, sx) == guide)){
						double dist = distance(Target, Source, tx, ty, sx, sy, PatchSize, bestDistance);

						if (bestDistance > dist){
							Offset.at<Vec2f>(ty, tx) = Vec2f(sx - tx, sy - ty);
							bestDistance = dist;
						}
					}
				}
				else if (!Hole.at<uchar>(sy, sx)){
					if (guide == 0 || Guide.at<uchar>(sy, sx) == guide){
						double dist = distance(Target, Source, tx, ty, sx, sy, PatchSize, bestDistance);
						if (dist >= 0){
							if (bestDistance > dist){
								Offset.at<Vec2f>(ty, tx) = offset;
								bestDistance = dist;
							}
						}
					}
				}
			}
			if (ty != yfirst){
				offset = Offset.at<Vec2f>(ty - step, tx);
				sx = tx + offset[0], sy = ty + offset[1];
				if (sy<0 || sy>Hole.rows - 1){
					sy = sy - step;
					if (sx<0 || sx>Hole.cols - 1 || sy<0 || sy>Hole.rows - 1){
						fprintf(stderr, "Wrong!%d,%d->%d,%d\n", tx, ty, sx, sy);
					}

					if (!Hole.at<uchar>(sy, sx) && (guide == 0 || Guide.at<uchar>(sy, sx) == guide)){
						double dist = distance(Target, Source, tx, ty, sx, sy, PatchSize, bestDistance);
						if (bestDistance > dist){
							Offset.at<Vec2f>(ty, tx) = Vec2f(sx - tx, sy - ty);
							bestDistance = dist;
						}
					}
				}
				else if (!Hole.at<uchar>(sy, sx)){
					if (guide == 0.1 || Guide.at<uchar>(sy, sx) == guide){
						double dist = distance(Target, Source, tx, ty, sx, sy, PatchSize, bestDistance);
						if (dist > 0){
							if (bestDistance > dist){
								Offset.at<Vec2f>(ty, tx) = offset;
								bestDistance = dist;
							}
						}
					}
				}
			}
			//random search
			for (float r = searchWindow; r >= 1; r *= alpha){
				int minx = tx + Offset.at<Vec2f>(ty, tx)[0] - r;
				int maxx = tx + Offset.at<Vec2f>(ty, tx)[0] + r;
				int miny = ty + Offset.at<Vec2f>(ty, tx)[1] - r;
				int maxy = ty + Offset.at<Vec2f>(ty, tx)[1] + r;
				if (searchSpaces[guide].size() == 0){
					guide = 0;
				}
				Vec2f s = searchSpaces[guide][rand() % searchSpaces[guide].size()];
				sx = s[0], sy = s[1];
				if (sx > maxx || sx<minx || sy>maxy || sy < miny){
					continue;
				}
				double dist = distance(Target, Source, tx, ty, sx, sy, PatchSize, bestDistance);
				if (dist >= 0){
					if (bestDistance > dist){
						bestDistance = dist;
						Offset.at<Vec2f>(ty, tx) = s - Vec2f(tx, ty);
					}
				}
			}
			Dist.at<double>(ty, tx) = bestDistance;
			Vec2f checks = Offset.at<Vec2f>(ty, tx) + Vec2f(tx, ty);
			if (checks[0] < 0 || checks[0] >= Source.cols || checks[1] < 0 || checks[1] >= Source.rows){
				fprintf(stdout, "not valid here (%d,%d)->(%f,%f):(%d,%d)", tx, ty, checks[0], checks[1], sx, sy);
			}

		}
	}
	return true;
}
void retarget(Mat SourceImage, Mat& TargetImage, Mat Hole, Mat Guide, Mat Const, Mat Offset, int minSize, int PatchSize, int searchWindowFacktor, int EMiterations, int ANNIterations){
	srand((unsigned int) (time(0)));
	int t1 = clock();
	Mat SrcImg;
	SourceImage.convertTo(SrcImg, CV_32FC3, 1 / 255.f);
	Mat DstImg;
	TargetImage.convertTo(DstImg, CV_32FC3, 1 / 255.f);
	//prepare space for pyramid
	int pyrHeight = 1;
	while (SrcImg.cols > minSize && SrcImg.rows > minSize){
		minSize *= 2;
		pyrHeight += 1;
	}
	Mat* Srcs = new Mat[pyrHeight];
	Mat* Dsts = new Mat[pyrHeight];
	Mat* Holes = new Mat[pyrHeight];
	Mat* Guides = new Mat[pyrHeight];
	Mat* Offs = new Mat[pyrHeight];
	Mat* Consts = new Mat[pyrHeight];
	Mat* Dists = new Mat[pyrHeight];
	Srcs[0] = SrcImg;
	Dsts[0] = DstImg;
	Holes[0] = Hole;
	Guides[0] = Guide;
	Offs[0] = Offset;
	Consts[0] = Const;
	Dists[0] = Mat(Offset.rows, Offset.cols, CV_64FC1, Scalar::all(-1));
	//pyrDown
	for (int i = 0; i < pyrHeight - 1; i++){
		pyrDown(Srcs[i], Srcs[i + 1]);
		pyrDown(Dsts[i], Dsts[i + 1]);
		pyrDownHole(Holes[i], Guides[i], Holes[i + 1], Guides[i + 1]);
		pyrDownOffs(Offs[i], Consts[i], Offs[i + 1], Consts[i + 1]);
		Dists[i + 1] = Mat(Offs[i + 1].rows, Offs[i + 1].cols, CV_64FC1, Scalar::all(-1));
		showAll(SrcImg, Srcs[i], Dsts[i], Holes[i], Guides[i], Offs[i], Consts[i]);
	}

	//release some levels of constraints
	for (int i = 0; i < pyrHeight / 2; i++){
		Consts[i] = Mat(Consts[i].rows, Consts[i].cols, Consts[i].type(), Scalar::all(0));
	}

	vector<vector<Vec2f>> searchSpaces;
	searchSpaces.reserve(16);
	searchSpaces.resize(1);
	//get search spaces for level[pyrHeight-1]
	for (int i = 0; i < Guides[pyrHeight - 1].rows; i++){
		for (int j = 0; j < Guides[pyrHeight - 1].cols; j++){
			if (Holes[pyrHeight - 1].at<uchar>(i, j)){
				continue;
			}
			int guide = Guides[pyrHeight - 1].at<uchar>(i, j);
			if (searchSpaces.size() < guide + 1){
				searchSpaces.resize(guide + 1);
			}
			if (guide != 0){
				searchSpaces[guide].push_back(Vec2f(j, i));
			}
			searchSpaces[0].push_back(Vec2f(j, i));
		}
	}
	//random init offset
	for (int i = 0; i < Offs[pyrHeight - 1].rows; i++){
		for (int j = 0; j < Offs[pyrHeight - 1].cols; j++){
			if (Consts[pyrHeight - 1].at<uchar>(i, j)){
				continue;
			}
			if (Holes[pyrHeight - 1].at<uchar>(i, j) == false){
				Vec2f s = searchSpaces[0][rand() % searchSpaces[0].size()];
				Offs[pyrHeight - 1].at<Vec2f>(i, j) = s - Vec2f(j, i);
			}
			else{
				int guide = Guides[pyrHeight - 1].at<uchar>(i, j);
				if (guide >= searchSpaces.size() || searchSpaces[guide].size() == 0){
					guide = 0;
				}
				Vec2f s = searchSpaces[guide][rand() % searchSpaces[guide].size()];
				Offs[pyrHeight - 1].at<Vec2f>(i, j) = s - Vec2f(j, i);
			}
		}
	}

	//pyrUp
	for (int i = pyrHeight - 1; i > 0; i--){
		fprintf(stdout, "pyrUp step %d\n", i);
		showAll(SrcImg, Srcs[i], Dsts[i], Holes[i], Guides[i], Offs[i], Consts[i]);
#ifdef _PATCHTABLE_
		PatchTable<float, float>* table = build_table(Srcs[i], Holes[i]);
#endif
		for (int j = 0; j < EMiterations; j++){
			fprintf(stdout, "EMiteration %d\n", j);
#ifndef _PATCHTABLE_
			//Ann iterations
			for (int k = 0; k < ANNIterations; k++){
				fprintf(stdout, "ANNiteration %d\n", k);
				Ann(Dsts[i], Srcs[i], Offs[i], Dists[i], Holes[i], Guides[i], Consts[i], searchSpaces, PatchSize, i*searchWindowFacktor, 0.5, k);
				showAll(SrcImg, Srcs[i], Dsts[i], Holes[i], Guides[i], Offs[i], Consts[i]);
			}
			Dists[i] = Mat(Dists[i].rows, Dists[i].cols, CV_64FC1, Scalar::all(-1));
#endif
			//patchtable
#ifdef _PATCHTABLE_
			table_match(table, Dsts[i], Srcs[i], Offs[i], Holes[i], PatchSize);
#endif
			generateImage(Offs[i], Srcs[i], Dsts[i], PatchSize);
			showAll(SrcImg, Srcs[i], Dsts[i], Holes[i], Guides[i], Offs[i], Consts[i]);
			showGUI(DstImg, Dsts[i]);
		}
#ifdef _PATCHTABLE_
		delete table;
#endif
		showGUI(DstImg, Dsts[i]);
		//get search space for next level
		searchSpaces.clear();
		searchSpaces.resize(1);
		for (int m = 0; m < Guides[i - 1].rows; m++){
			for (int n = 0; n < Guides[i - 1].cols; n++){
				if (Holes[i - 1].at<uchar>(m, n)){
					continue;
				}
				int guide = Guides[i - 1].at<uchar>(m, n);
				if (searchSpaces.size() < guide + 1){
					searchSpaces.resize(guide + 1);
				}
				if (guide != 0){
					searchSpaces[guide].push_back(Vec2f(n, m));
				}
				searchSpaces[0].push_back(Vec2f(n, m));
			}
		}
		//pryUp offset and generator target image for next level
		pyrUpOffset(Offs[i], Offs[i - 1], Holes[i - 1], Guides[i - 1], Consts[i - 1], searchSpaces);
		showAll(SrcImg, Srcs[i], Dsts[i], Holes[i], Guides[i], Offs[i], Consts[i]);
		generateImage(Offs[i - 1], Srcs[i - 1], Dsts[i - 1], PatchSize, true);
		showAll(SrcImg, Srcs[i], Dsts[i], Holes[i], Guides[i], Offs[i], Consts[i]);
	}
	fprintf(stdout, "pyrUp step 0\n");
	//generate final DstImg
	DstImg = Dsts[0];
	PatchTable<float, float>* table = build_table(Srcs[0], Holes[0]);
	for (int j = 0; j < EMiterations; j++){
		fprintf(stdout, "EMiteration %d\n", j);
		//Ann iterations
		/*for (int k = 0; k < ANNIterations; k++){
			fprintf(stdout, "ANNiteration %d\n", k);
			Ann(DstImg, Srcs[0], Offs[0], Dists[0], Holes[0], Guides[0], Consts[0], searchSpaces, PatchSize, 1, 0.5, k);
			}
			Dists[0] = Mat(Dists[0].rows, Dists[0].cols, CV_64FC1, Scalar::all(-1));*/
		table_match(table, Dsts[0], Srcs[0], Offs[0], Holes[0], PatchSize);
		generateImage(Offs[0], Srcs[0], DstImg, PatchSize);
	}
	delete table;
	delete[] Srcs;
	delete[] Dsts;
	delete[] Holes;
	delete[] Guides;
	delete[] Dists;
	//output DstImg as Dst.png & Offset map as Ann.png
	imwrite("Dst.png", TargetImage);
	int row = Offs[0].rows, col = Offs[0].cols;
	Mat OffImg(row, col, CV_8UC3);
	for (int i = 0; i < row; i++){
		for (int j = 0; j < col; j++){
			int r = (Offs[0].at<Vec2f>(i, j)[0] + (SrcImg.cols + DstImg.cols)) * 255 / (2 * (SrcImg.cols + DstImg.cols));
			int g = (Offs[0].at<Vec2f>(i, j)[1] + (SrcImg.rows + DstImg.rows)) * 255 / (2 * (SrcImg.rows + DstImg.rows));
			OffImg.at<Vec3b>(i, j) = Vec3b(0, g, r);
		}
	}
	imwrite("Ann.png", OffImg);
	delete[] Offs;
	DstImg.convertTo(TargetImage, CV_8UC3, 255.f);
	int t2 = clock();
	fprintf(stdout, "Cost %d ms.\n", t2 - t1);
}

double distance(Mat Target, Mat Source, int tx, int ty, int sx, int sy, int PatchSize, double abort){
	PatchSize /= 2;
	int left = MIN(tx, PatchSize), right = MIN(Target.cols - tx, PatchSize + 1);
	int up = MIN(ty, PatchSize), down = MIN(Target.rows - ty, PatchSize + 1);
	if (sx - left<0 || sx + right>Source.cols || sy - up<0 || sy + down>Source.rows){
		return -1;
	}
	double dist = 0;
	abort *= (right + left - 1)*(down + up - 1) * 3;
	for (int i = -left; i < right; i++){
		for (int j = -up; j < down; j++){
			Vec3f d = Target.at<Vec3f>(ty + j, tx + i) - Source.at<Vec3f>(sy + j, sx + i);
			dist += d.dot(d);
		}
		if (dist > abort){
			return MAX_DIST;
		}
	}
	dist = dist / (right + left - 1) / (down + up - 1) / 3;
	return dist;
}
void pyrDownHole(Mat HoleIn, Mat GuideIn, Mat& HoleOut, Mat& GuideOut){
	int rows = (HoleIn.rows + 1) / 2;
	int cols = (HoleIn.cols + 1) / 2;
	HoleOut = Mat(rows, cols, HoleIn.type());
	GuideOut = Mat(rows, cols, GuideIn.type());
	for (int i = 0; i < rows; i++){
		for (int j = 0; j < cols; j++){
			if (HoleIn.at<uchar>(2 * i, 2 * j)){
				HoleOut.at<uchar>(i, j) = HoleIn.at<uchar>(2 * i, 2 * j);
				GuideOut.at<uchar>(i, j) = GuideIn.at<uchar>(2 * i, 2 * j);
			}
			else if (HoleIn.rows>2 * i + 1 && HoleIn.at<uchar>(2 * i + 1, 2 * j)){
				HoleOut.at<uchar>(i, j) = HoleIn.at<uchar>(2 * i + 1, 2 * j);
				GuideOut.at<uchar>(i, j) = GuideIn.at<uchar>(2 * i + 1, 2 * j);
			}
			else if (HoleIn.cols>2 * j + 1 && HoleIn.at<uchar>(2 * i, 2 * j + 1)){
				HoleOut.at<uchar>(i, j) = HoleIn.at<uchar>(2 * i, 2 * j + 1);
				GuideOut.at<uchar>(i, j) = GuideIn.at<uchar>(2 * i, 2 * j + 1);
			}
			else if (HoleIn.rows > 2 * i + 1 && HoleIn.cols > 2 * j + 1 && HoleIn.at<uchar>(2 * i + 1, 2 * j + 1)){
				HoleOut.at<uchar>(i, j) = HoleIn.at<uchar>(2 * i + 1, 2 * j + 1);
				GuideOut.at<uchar>(i, j) = GuideIn.at<uchar>(2 * i + 1, 2 * j + 1);
			}
			else{
				HoleOut.at<uchar>(i, j) = 0;
				GuideOut.at<uchar>(i, j) = GuideIn.at<uchar>(2 * i, 2 * j);
			}
		}
	}
}
void pyrDownOffs(Mat OffIn, Mat ConstIn, Mat& OffOut, Mat& ConstOut){
	int rows = (OffIn.rows + 1) / 2;
	int cols = (OffIn.cols + 1) / 2;
	OffOut = Mat(rows, cols, OffIn.type());
	ConstOut = Mat(rows, cols, ConstIn.type());
	for (int i = 0; i < rows; i++){
		for (int j = 0; j < cols; j++){
			if (ConstIn.at<uchar>(2 * i, 2 * j)){
				ConstOut.at<uchar>(i, j) = ConstIn.at<uchar>(2 * i, 2 * j);
				OffOut.at<Vec2f>(i, j) = OffIn.at<Vec2f>(2 * i, 2 * j) / 2;
			}
			else if (ConstIn.rows>2 * i + 1 && ConstIn.at<uchar>(2 * i + 1, 2 * j)){
				ConstOut.at<uchar>(i, j) = ConstIn.at<uchar>(2 * i + 1, 2 * j);
				OffOut.at<Vec2f>(i, j) = OffIn.at<Vec2f>(2 * i + 1, 2 * j) / 2;
			}
			else if (ConstIn.cols>2 * j + 1 && ConstIn.at<uchar>(2 * i, 2 * j + 1)){
				ConstOut.at<uchar>(i, j) = ConstIn.at<uchar>(2 * i, 2 * j + 1);
				OffOut.at<Vec2f>(i, j) = OffIn.at<Vec2f>(2 * i, 2 * j + 1) / 2;
			}
			else if (ConstIn.rows > 2 * i + 1 && ConstIn.cols > 2 * j + 1 && ConstIn.at<uchar>(2 * i + 1, 2 * j + 1)){
				ConstOut.at<uchar>(i, j) = ConstIn.at<uchar>(2 * i + 1, 2 * j + 1);
				OffOut.at<Vec2f>(i, j) = OffIn.at<Vec2f>(2 * i + 1, 2 * j + 1) / 2;
			}
			else{
				ConstOut.at<uchar>(i, j) = 0;
				OffOut.at<Vec2f>(i, j) = OffIn.at<Vec2f>(2 * i, 2 * j) / 2;
			}
		}
	}
}
void pyrUpOffset(Mat OffsetIn, Mat& OffsetOut, Mat HoleOut, Mat GuideOut, Mat ConstOut, vector<vector<Vec2f>> searchSpaces){
	for (int i = 0; i < OffsetOut.rows; i++){
		for (int j = 0; j < OffsetOut.cols; j++){
			if (ConstOut.at<uchar>(i, j)){
				continue;
			}
			OffsetOut.at<Vec2f>(i, j) = OffsetIn.at<Vec2f>(i / 2, j / 2) * 2;
			if (HoleOut.at<uchar>(i, j)){
				int guide = GuideOut.at<uchar>(i, j);
				if (guide >= searchSpaces.size() || searchSpaces[guide].size() == 0){
					guide = 0;
				}
				Vec2f offset = OffsetOut.at<Vec2f>(i, j);
				if (i + offset[1]<0 || i + offset[1]>GuideOut.rows - 1 || j + offset[0]<0 || j + offset[0]>GuideOut.cols - 1 ||
					GuideOut.at<uchar>(i + offset[1], j + offset[0]) != guide){
					OffsetOut.at<Vec2f>(i, j) = searchSpaces[guide][rand() % searchSpaces[guide].size()];
				}
			}
		}
	}
}
void generateImage(Mat Offset, Mat SrcImg, Mat& DstImg, int patchSize, bool up){
	patchSize = patchSize / 2;
	Mat Count(DstImg.rows, DstImg.cols, CV_32FC1, Scalar::all(0));
	DstImg = Mat(DstImg.rows, DstImg.cols, CV_32FC3, Scalar::all(0));
	int fix = 1;
	if (up){
		patchSize *= 2;
		fix = 2;
	}
	for (int tx = patchSize; tx < DstImg.cols - patchSize - fix; tx++){
		for (int ty = patchSize; ty < DstImg.rows - patchSize - fix; ty++){
			Vec2f offset = Offset.at<Vec2f>(ty, tx);
			int sx = offset[0] + tx, sy = offset[1] + ty;
			int left = MIN(tx, patchSize), right = MIN(DstImg.cols - tx, patchSize + fix + 1);
			int up = MIN(ty, patchSize), down = MIN(DstImg.rows - ty, patchSize + fix + 1);
			if (sx - left<0 || sx + right>SrcImg.cols || sx - left > sx + right || sy - up<0 || sy + down>SrcImg.rows || sy - up > sy + down){
//				fprintf(stderr, "Warning:target(%d,%d) match source(%d,%d) is not valid.\n", tx, ty, sx, sy);
			}
			else{
				Mat Patch = SrcImg(Range(sy - up, sy + down), Range(sx - left, sx + right));
				Mat Dst = DstImg(Range(ty - up, ty + down), Range(tx - left, tx + right));
				Mat Cou = Count(Range(ty - up, ty + down), Range(tx - left, tx + right));
				Dst += Patch;
				Cou += 1;
			}
		}
	}
	std::vector<Mat> channels(3);
	split(DstImg, channels);
	channels[0] /= Count;
	channels[1] /= Count;
	channels[2] /= Count;
	merge(channels, DstImg);
}
