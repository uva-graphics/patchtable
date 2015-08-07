#include "PatchMatch.h"
#include "patchtable.h"
#include<ctime>
#define MAX_DIST 1000000 //maximum distance for the distance function
#define ZERO 0.00000000001
//#define _DEBUG_
//#define _FIX_
//for debug
bool final_flag = false;
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
void pyrDownCover(Mat CoverIn, Mat& CoverOut);
void pyrUpOffset(Mat OffsetIn, Mat& OffsetOut, Mat HoleOut, Mat GuideOut, Mat ConstOut, vector<vector<Vec2f>> searchSpaces);
void generateImage(Mat Offset, Mat SrcImg, Mat& DstImg, Mat Const, int patchSize, bool up = false);
void getCombination(Mat& combine, Mat& hole, Mat& guide);
void getTransformation(Mat& transformation, Mat OriImg, Mat& hole, Mat Cover, Mat& guide);
void refineImage(Mat& Src, Mat& Dst, Mat& Const, int PatchSize);
void assessResult(Mat& Offset, Mat& Const, float& coherence);

PatchTable<float, float>* build_table(Mat B, Mat Hole, int coherence_spatial, int coherence_temporal){
	PatchTableParams *p = new PatchTableParams();
	//p->set_speed(0);
	p->coherence_spatial = coherence_spatial;
	p->coherence_temporal = coherence_temporal;
	p->calc_exact_dist = true;
	p->limit = 1000000;
	Mat B_;
	Array<int32_t> allow(Hole.rows, Hole.cols);
	//fprintf(stdout, "%d %d %d %d\n", B.rows,B.cols,Hole.rows,Hole.cols);
	for (int row = 0; row < Hole.rows; row++){
		for (int col = 0; col < Hole.cols; col++){
			if (Hole.at<uchar>(row, col) == 1){
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

void simple_init(Mat& Dst, Mat Src, Mat Cover, Mat Hole, int PatchSize)
{
	//get the points to be covered
	vector<Point> cover_pnts;
	for (int row = 0; row < Cover.rows; ++row) {
		for (int col = 0; col < Cover.cols; ++col) {
			if (Cover.at<uchar>(row, col) == 1) {
				cover_pnts.push_back(Point(row, col));
			}
		}
	}
	//for every points to be covered
	for (size_t i = 0; i < cover_pnts.size(); ++i) {
		int x = cover_pnts[i].x, y = cover_pnts[i].y;
		if (x < PatchSize - 1 || y < PatchSize - 1) {
			return;
		}
		//get the color of the patch
		vector<Vec3f> cover_rgbs;
		for (int row = x; row > x - PatchSize; --row) {
			for (int col = y; col > y - PatchSize; --col) {
				if (x != row || y != col) {
					cover_rgbs.push_back(Dst.at<Vec3f>(row, col));
				}
			}
		}
		//search for the best patch
		float minDis = 1e5;
		Point BestPatch(PatchSize - 1, PatchSize - 1);
		for (int px = PatchSize - 1; px < Src.rows; ++px) {
			for (int py = PatchSize - 1; py < Src.cols; ++py) {
				bool hole_flag = false;
				for (int row = px; row > px - PatchSize; --row) {
					for (int col = py; col > py - PatchSize; --col) {
						if (!Hole.at<uchar>(row, col)) {
							hole_flag = true;
							break;
						}
					}
					if (hole_flag) {
						break;
					}
				}
				if (hole_flag) {
					continue;
				}
				float dis = 0;
				size_t iter = 0;
				for (int row = px; row > px - PatchSize; --row) {
					for (int col = py; col > py - PatchSize; --col) {
						if (px != row || py != col) {
							Vec3f vec = cover_rgbs[iter] - Src.at<Vec3f>(row, col);
							dis += vec.dot(vec);
							++iter;
						}
					}
				}
				if (dis < minDis) {
					BestPatch.x = px; BestPatch.y = py;
					minDis = dis;
				}
			}
		}
		Dst.at<Vec3f>(x, y) = Src.at<Vec3f>(BestPatch.x, BestPatch.y);
		//fprintf(stdout, "%f  (%d, %d)\n", minDis, BestPatch.x, BestPatch.y);
	}
}

bool isBoundry(Mat& Cover, int row, int col)
{
	if (!Cover.at<uchar>(row, col)) {
		if (row + 1 < Cover.rows && Cover.at<uchar>(row + 1, col)) {
			return true;
		}
		if (row - 1 >= 0 && Cover.at<uchar>(row - 1, col)) {
			return true;
		}
		if (col + 1 < Cover.cols && Cover.at<uchar>(row, col + 1)) {
			return true;
		}
		if (col - 1 >= 0 && Cover.at<uchar>(row, col - 1)) {
			return true;
		}
	}
	return false;
}

void interpolate_init(Mat& Dst, Mat& Src, Mat& Cover, Mat& Hole, int PatchSize)
{
	//int ratio = 16;
	//Mat _cover(Cover.rows*ratio, Cover.cols*ratio, Cover.type());
	//for (int row = 0; row < Cover.rows; ++row) {
	//	for (int col = 0; col < Cover.cols; ++col) {
	//		for (int i = 0; i < ratio; ++i) {
	//			for (int j = 0; j < ratio; ++j) {
	//				_cover.at<uchar>(row*ratio + i, col*ratio + j) = Cover.at<uchar>(row, col) * 255;
	//			}
	//		}
	//	}
	//}
	//imshow("cover", _cover);
	//waitKey(0);
	//get the points to be covered
	vector<Point> cover_pnts;
	for (int row = 0; row < Cover.rows; ++row) {
		for (int col = 0; col < Cover.cols; ++col) {
			if (Cover.at<uchar>(row, col) == 1) {
				cover_pnts.push_back(Point(row, col));
			}
		}
	}
	if (cover_pnts.empty()) {
		return;
	}
	//get the boundry
	vector<Point> boundry;
	bool start = false;
	for (int row = 0; row < Cover.rows; ++row) {
		for (int col = 0; col < Cover.cols; ++col) {
			if (isBoundry(Cover, row, col)) {
				boundry.push_back(Point(row, col));
				start = true;
				break;
			}
		}
		if (start) {
			break;
		}
	}
	int direction = 0;
	while (true) {
		int row = boundry[boundry.size() - 1].x, col = boundry[boundry.size() - 1].y;
		if (boundry.size() > 1 && row == boundry[0].x && col == boundry[0].y) {
			boundry.pop_back();
			break;
		}
		int _row = row, _col = col;
		if (boundry.size() > 1) {
			_row = boundry[boundry.size() - 2].x;
			_col = boundry[boundry.size() - 2].y;
		}
		bool flag = false;
		for (int x = MIN(row + 1, Cover.rows - 1); x >= MAX(row - 1, 0); --x) {
			for (int y = MIN(col + 1, Cover.cols - 1); y >= MAX(col - 1, 0); --y) {
				if ((x != row || y != col) && (x != _row || y != _col)) {
					if (isBoundry(Cover, x, y)) {
						boundry.push_back(Point(x, y));
						flag = true;
						break;
					}
				}
			}
			if (flag) {
				break;
			}
		}
	}
	//for every points to be covered
	for (size_t i = 0; i < cover_pnts.size(); ++i){
		vector<double> weights;
		double weight_sum = 0.0;
		for (size_t j = 0; j < boundry.size(); ++j) {
			Vec2d x(cover_pnts[i].x, cover_pnts[i].y);
			Vec2d p(boundry[j].x, boundry[j].y);
			Vec2d p1, p2;
			if (j == 0) {
				p1 = Vec2d(boundry[boundry.size() - 1].x, boundry[boundry.size() - 1].y);
				p2 = Vec2d(boundry[j + 1].x, boundry[j + 1].y);
			}
			else if (j == boundry.size() - 1) {
				p1 = Vec2d(boundry[j - 1].x, boundry[j - 1].y);
				p2 = Vec2d(boundry[0].x, boundry[j - 1].y);
			}
			else {
				p1 = Vec2d(boundry[j - 1].x, boundry[j - 1].y);
				p2 = Vec2d(boundry[j + 1].x, boundry[j + 1].y);
			}
			Vec2d v1 = p1 - x, v = p - x, v2 = p2 - x;
			double cos1 = v1.dot(v) / sqrt(v1.dot(v1)*v.dot(v));
			double cos2 = v2.dot(v) / sqrt(v2.dot(v2)*v.dot(v));
			double weight = (sqrt((1 - cos1) / (1 + cos1)) + sqrt((1 - cos2) / (1 + cos2))) / sqrt(v.dot(v));
			weights.push_back(weight);
			weight_sum += weight;
		}
		Vec3f color(0, 0, 0);
		for (size_t j = 0; j < boundry.size(); ++j) {
			int row = boundry[j].x, col = boundry[j].y;
			color += (weights[j] / weight_sum) * Dst.at<Vec3f>(row, col);
		}
		if (color.val[0] > 1.0) {
			color.val[0] = 1;
		}
		if (color.val[1] > 1.0) {
			color.val[1] = 1;
		}
		if (color.val[2] > 1.0) {
			color.val[2] = 1;
		}
		Dst.at<Vec3f>(cover_pnts[i].x, cover_pnts[i].y) = color;
	}
	//int ratio = 16;
	//Mat _dst(Dst.rows*ratio, Dst.cols*ratio, Dst.type());
	//for (int row = 0; row < Dst.rows; ++row) {
	//	for (int col = 0; col < Dst.cols; ++col) {
	//		for (int i = 0; i < ratio; ++i) {
	//			for (int j = 0; j < ratio; ++j) {
	//				_dst.at<Vec3f>(row*ratio + i, col*ratio + j) = Dst.at<Vec3f>(row, col);
	//			}
	//		}
	//	}
	//}
	//imshow("dst", _dst);
	//waitKey(0);
}

void match_init(Mat& Dst, Mat& Src, Mat& Cover, Mat& Hole, Mat& Offset, int PatchSize, int matchIndex)
{
	int minRow = Cover.rows, minCol = Cover.cols, maxRow = 0, maxCol = 0;
	vector<Point> boundry;
	vector<Point> cover_pnts;
	for (int i = 0; i < Cover.rows; ++i) {
		for (int j = 0; j < Cover.cols; ++j) {
			if (!Cover.at<uchar>(i, j)) {
				if ((i - 1 >= 0 && Cover.at<uchar>(i - 1, j)) || (i + 1 < Cover.rows && Cover.at<uchar>(i + 1, j))
					|| (j - 1 >= 0 && Cover.at<uchar>(i, j - 1)) || (j + 1 < Cover.cols && Cover.at<uchar>(i, j + 1))) {
					boundry.push_back(Point(i, j));
					if (i < minRow) {
						minRow = i;
					}
					if (i > maxRow) {
						maxRow = i;
					}
					if (j < minCol) {
						minCol = j;
					}
					if (j > maxCol) {
						maxCol = j;
					}
				}
			}
			else {
				cover_pnts.push_back(Point(i, j));
			}
		}
	}

	vector<float> distances;
	vector<int> rowOffs, colOffs;
	for (int row = 0; row < Src.rows - maxRow + minRow; ++row) {
		for (int col = 0; col < Src.cols - maxCol + minCol; ++col) {
			int rowOff = row - minRow, colOff = col - minCol;
			float dis = 0.0;
			bool holeFlag = false;
			for (size_t i = 0; i < cover_pnts.size(); ++i) {
				if (!Hole.at<uchar>(cover_pnts[i].x + rowOff, cover_pnts[i].y + colOff)) {
					holeFlag = true;
					break;
				}
			}
			if (holeFlag) {
				continue;
			}
			for (size_t i = 0; i < boundry.size(); ++i) {
				int x = boundry[i].x, y = boundry[i].y;
				Vec3f diff = Dst.at<Vec3f>(x, y) - Src.at<Vec3f>(x + rowOff, y + colOff);
				dis += diff.dot(diff);
			}
			if (holeFlag) {
				continue;
			}
			if (distances.size() < 5) {
				distances.push_back(dis);
				rowOffs.push_back(rowOff);
				colOffs.push_back(colOff);
			}
			else {
				float maxDis = 0;
				size_t idx = 0;
				for (size_t i = 0; i < 5; ++i) {
					if (distances[i] > maxDis) {
						maxDis = distances[i];
						idx = i;
					}
				}
				if (maxDis > dis) {
					distances[idx] = dis;
					rowOffs[idx] = rowOff;
					colOffs[idx] = colOff;
				}
			}
		}
	}

	for (int i = 0; i < Dst.rows; ++i) {
		for (int j = 0; j < Dst.cols; ++j) {
			if (Cover.at<uchar>(i, j)) {
				Dst.at<Vec3f>(i, j) = Src.at<Vec3f>(i + rowOffs[matchIndex], j + colOffs[matchIndex]);
				Offset.at<Vec2f>(i, j) = Vec2f(colOffs[matchIndex], rowOffs[matchIndex]);
			}
		}
	}

}
void match_init1(Mat& Dst, Mat& Src, Mat& Cover, Mat& Hole, Mat& Offset, int PatchSize, int matchIndex)
{
	int leftCol = 0, rightCol = 0;
	bool coverFlag = false;
	for (int col = 0; col < Cover.cols; ++col) {
		if (!coverFlag && Cover.at<uchar>(0, col)) {
			leftCol = col;
			coverFlag = true;
		}
		if (coverFlag && !Cover.at<uchar>(0, col)) {
			rightCol = col - 1;
			break;
		}
	}

	vector<float> distances;
	vector<int> colOffs;
	for (int col = 0; col < Src.cols - rightCol + leftCol; ++col) {
		int colOff = col - leftCol;
		float dis1 = 0.0, dis2 = 0.0;
		for (int row = 0; row < Src.rows; ++row) {
			Vec3f diff1 = Dst.at<Vec3f>(row, leftCol - 1) - Src.at<Vec3f>(row, col);
			Vec3f diff2 = Dst.at<Vec3f>(row, rightCol + 1) - Src.at<Vec3f>(row, col + rightCol - leftCol);
			dis1 += diff1.dot(diff1);
			dis2 += diff2.dot(diff2);
		}
		float dis = dis1 + dis2;
		if (distances.size() < 5) {
			distances.push_back(dis);
			colOffs.push_back(colOff);
		}
		else {
			float maxDis = 0;
			size_t idx = 0;
			for (size_t i = 0; i < 5; ++i) {
				if (distances[i] > maxDis) {
					maxDis = distances[i];
					idx = i;
				}
			}
			if (maxDis > dis) {
				distances[idx] = dis;
				colOffs[idx] = colOff;
			}
		}
	}

	for (int i = 0; i < Dst.rows; ++i) {
		for (int j = 0; j < Dst.cols; ++j) {
			if (Cover.at<uchar>(i, j)) {
				Dst.at<Vec3f>(i, j) = Src.at<Vec3f>(i, j + colOffs[matchIndex]);
				Offset.at<Vec2f>(i, j) = Vec2f(colOffs[matchIndex], 0);
			}
		}
	}

}

void table_match(PatchTable<float, float> *table, Mat A, Mat B, Mat& offset, Mat Hole, int PatchSize){
	double ratio = 0.8;
	Array<double> ann;
	Mat A_;
	A.convertTo(A_, CV_8UC3, 255.f);
	imwrite("A.png", A_);
	Array<float> a(load_color_image<float>("A.png"));
	double latest_time = 1e100;
//	for (int iter = 0; iter < 5; iter++){               // Connelly: replaced 5 iterations of lookup() with one iteration.
    printf("Before lookup()\n");
		latest_time = table->lookup(a, ann,&ann);
//	}
    printf("After lookup()\n");
	for (int i = 0; i < offset.rows - 7; i++){
		for (int j = 0; j < offset.cols - 7; j++){
			int tx = j + 3, ty = i + 3;
			Vec2f off = offset.at<Vec2f>(ty, tx);
			int sx = tx + off[0], sy = ty + off[1];
			double ori_dist = distance(A, B, tx, ty, sx, sy, PatchSize);
			float x = float(ann(i, j, NNF_X));
			float y = float(ann(i, j, NNF_Y));
			double dist = distance(A, B, tx, ty, x + 3, y + 3, PatchSize);
			if (sy < Hole.rows && sx < Hole.cols && Hole.at<uchar>(sy, sx) && dist < ori_dist*ratio){
				offset.at<Vec2f>(ty, tx) = Vec2f(x - j, y - i);
			}
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
void retarget(Mat SourceImage, Mat& TargetImage, Mat Cover, Mat Hole, Mat Guide, Mat Const, Mat Offset, 
	int minSize, int PatchSize, int searchWindowFacktor, int EMiterations, int ANNIterations, int coherence_spatial, int coherence_temporal, int matchIndex, float& coherence){
	srand((unsigned int) (time(0)));
	int t1 = clock();
	Mat SrcImg;
	SourceImage.convertTo(SrcImg, CV_32FC3, 1 / 255.f);
	Mat DstImg;
	TargetImage.convertTo(DstImg, CV_32FC3, 1 / 255.f);
	//prepare space for pyramid
	int pyrHeight = 1;
	while (DstImg.cols > minSize && DstImg.rows > minSize){
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
	Mat* Covers = new Mat[pyrHeight];
	getCombination(Srcs[0], Holes[0], Guides[0]);
	//getTransformation(Srcs[0], TargetImage, Holes[0], Cover, Guides[0]);
	//Srcs[0] = SrcImg;
	Dsts[0] = DstImg;
	Mat Ori_Dst(DstImg.rows, DstImg.cols, DstImg.type());
	DstImg.copyTo(Ori_Dst);
	//Holes[0] = Hole;
	//Guides[0] = Guide;
	Offs[0] = Offset;
	Consts[0] = Const;
	Covers[0] = Cover;
	Dists[0] = Mat(Offset.rows, Offset.cols, CV_64FC1, Scalar::all(-1));
	//pyrDown
	for (int i = 0; i < pyrHeight - 1; i++){
		pyrDown(Srcs[i], Srcs[i + 1]);
		pyrDown(Dsts[i], Dsts[i + 1]);
		pyrDownCover(Covers[i], Covers[i + 1]);
		pyrDownHole(Holes[i], Guides[i], Holes[i + 1], Guides[i + 1]);
		pyrDownOffs(Offs[i], Consts[i], Offs[i + 1], Consts[i + 1]);
		Dists[i + 1] = Mat(Offs[i + 1].rows, Offs[i + 1].cols, CV_64FC1, Scalar::all(-1));
		showAll(SrcImg, Srcs[i], Dsts[i], Holes[i], Guides[i], Offs[i], Consts[i]);
	}

	//release some levels of constraints
	//for (int i = 0; i < pyrHeight / 2; i++){
	//	Consts[i] = Mat(Consts[i].rows, Consts[i].cols, Consts[i].type(), Scalar::all(0));
	//}

	vector<vector<Vec2f>> searchSpaces;
	searchSpaces.reserve(16);
	searchSpaces.resize(1);
	//get search spaces for level[pyrHeight-1]
	for (int i = 0; i < Guides[pyrHeight - 1].rows; i++){
		for (int j = 0; j < Guides[pyrHeight - 1].cols; j++){
			if (!Holes[pyrHeight - 1].at<uchar>(i, j)){
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
			//if (Consts[pyrHeight - 1].at<uchar>(i, j)){
			//	continue;
			//}
			//if (Holes[pyrHeight - 1].at<uchar>(i, j) == false){
				Vec2f s = searchSpaces[0][rand() % searchSpaces[0].size()];
				Offs[pyrHeight - 1].at<Vec2f>(i, j) = s - Vec2f(j, i);
			//}
			//else{
			//	int guide = Guides[pyrHeight - 1].at<uchar>(i, j);
			//	if (guide >= searchSpaces.size() || searchSpaces[guide].size() == 0){
			//		guide = 0;
			//	}
			//	Vec2f s = searchSpaces[guide][rand() % searchSpaces[guide].size()];
			//	Offs[pyrHeight - 1].at<Vec2f>(i, j) = s - Vec2f(j, i);
			//}
		}
	}

	match_init(Dsts[pyrHeight - 1], Srcs[pyrHeight - 1], Covers[pyrHeight - 1], Holes[pyrHeight - 1], Offs[pyrHeight - 1], PatchSize, matchIndex);
	//interpolate_init(Dsts[pyrHeight - 1], Srcs[pyrHeight - 1], Covers[pyrHeight - 1], Holes[pyrHeight - 1], PatchSize);
	std::ofstream fout("time.txt");

	//pyrUp
	int iter = 0;
	for (int i = pyrHeight - 1; i > 0; i--){
		++iter;
		fprintf(stdout, "pyrUp step %d\n", i);
		showAll(SrcImg, Srcs[i], Dsts[i], Holes[i], Guides[i], Offs[i], Consts[i]);
		int time1 = clock();
		PatchTable<float, float>* table = build_table(Srcs[i], Holes[i], coherence_spatial, coherence_temporal);
		int time2 = clock();
		fout << "Pyramid iteration " << iter << ":" << std::endl;
		fout << "Building patchtable costs " << time2 - time1 << "ms." << std::endl;
		for (int j = 0; j < EMiterations; j++){
			fprintf(stdout, "EMiteration %d\n", j);
			//Ann iterations
			/*for (int k = 0; k < ANNIterations; k++){
				fprintf(stdout, "ANNiteration %d\n", k);
				Ann(Dsts[i], Srcs[i], Offs[i], Dists[i], Holes[i], Guides[i], Consts[i], searchSpaces, PatchSize, i*searchWindowFacktor, 0.5, k);
				showAll(SrcImg, Srcs[i], Dsts[i], Holes[i], Guides[i], Offs[i], Consts[i]);
				}
				Dists[i] = Mat(Dists[i].rows, Dists[i].cols, CV_64FC1, Scalar::all(-1));*/
			//patchtable
			table_match(table, Dsts[i], Srcs[i], Offs[i], Holes[i], PatchSize);
			generateImage(Offs[i], Srcs[i], Dsts[i], Consts[i], PatchSize);
			showAll(SrcImg, Srcs[i], Dsts[i], Holes[i], Guides[i], Offs[i], Consts[i]);
			showGUI(DstImg, Dsts[i]);
		}

		delete table;
		showGUI(DstImg, Dsts[i]);
		//get search space for next level
		searchSpaces.clear();
		searchSpaces.resize(1);
		for (int m = 0; m < Guides[i - 1].rows; m++){
			for (int n = 0; n < Guides[i - 1].cols; n++){
				if (!Holes[i - 1].at<uchar>(m, n)){
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
		generateImage(Offs[i - 1], Srcs[i - 1], Dsts[i - 1], Consts[i - 1], PatchSize, true);
		showAll(SrcImg, Srcs[i], Dsts[i], Holes[i], Guides[i], Offs[i], Consts[i]);
		int time3 = clock();
		fout << "The rest time is " << time3 - time2 << "ms." << std::endl << std::endl;
	}

	fprintf(stdout, "pyrUp step 0\n");
	++iter;
	//generate final DstImg
	DstImg = Dsts[0];
	int time1 = clock();
	PatchTable<float, float>* table = build_table(Srcs[0], Holes[0], coherence_spatial, coherence_temporal);
	int time2 = clock();
	fout << "Pyramid iteration " << iter << ":" << std::endl;
	fout << "Building patchtable costs " << time2 - time1 << "ms." << std::endl;
	for (int j = 0; j < EMiterations; j++){
		fprintf(stdout, "EMiteration %d\n", j);
		//Ann iterations
		/*for (int k = 0; k < ANNIterations; k++){
			fprintf(stdout, "ANNiteration %d\n", k);
			Ann(DstImg, Srcs[0], Offs[0], Dists[0], Holes[0], Guides[0], Consts[0], searchSpaces, PatchSize, 1, 0.5, k);
			}
			Dists[0] = Mat(Dists[0].rows, Dists[0].cols, CV_64FC1, Scalar::all(-1));*/
		table_match(table, Dsts[0], Srcs[0], Offs[0], Holes[0], PatchSize);
		generateImage(Offs[0], Srcs[0], DstImg, Consts[0], PatchSize);
	}
	int time3 = clock();
	fout << "The rest time is " << time3 - time2 << "ms." << std::endl;

	//refineImage(Ori_Dst, DstImg, Consts[0], PatchSize);
	assessResult(Offs[0], Consts[0], coherence);

	delete table;
	delete[] Srcs;
	delete[] Dsts;
	delete[] Holes;
	delete[] Guides;
	delete[] Dists;
	DstImg.convertTo(TargetImage, CV_8UC3, 255.f);
	string result_name = "result_";
	result_name.append(std::to_string(matchIndex));
	result_name.append("_");
	result_name.append(std::to_string(coherence));
	result_name.append(".png");
	imwrite(result_name.c_str(), TargetImage);
	Mat OffImg(Offs[0].rows, Offs[0].cols, CV_8UC3);
	int row = Offs[0].rows, col = Offs[0].cols;
	for (int i = 0; i<Offs[0].rows; i++){
		for (int j = 0; j<Offs[0].cols; j++){
			int r = (Offs[0].at<Vec2f>(i*Offs[0].rows / row, j*Offs[0].cols / col)[0] + (col)) * 255 / (2 * col);
			int g = (Offs[0].at<Vec2f>(i*Offs[0].rows / row, j*Offs[0].cols / col)[1] + (row)) * 255 / (2 * row);
			OffImg.at<Vec3b>(i, j) = Vec3b(0, g, r);
		}
	}
	imwrite("Ann.png", OffImg);
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
void pyrDownCover(Mat CoverIn, Mat& CoverOut) {
	pyrDown(CoverIn, CoverOut);
	for (int i = 0; i < CoverOut.rows; ++i) {
		CoverOut.at<uchar>(i, 0) = 0;
		CoverOut.at<uchar>(i, CoverOut.cols - 1) = 0;
	}
	for (int i = 0; i < CoverOut.cols; ++i) {
		CoverOut.at<uchar>(0, i) = 0;
		CoverOut.at<uchar>(CoverOut.rows - 1, i) = 0;
	}
}
void pyrUpOffset(Mat OffsetIn, Mat& OffsetOut, Mat HoleOut, Mat GuideOut, Mat ConstOut, vector<vector<Vec2f>> searchSpaces){
	for (int i = 0; i < OffsetOut.rows; i++){
		for (int j = 0; j < OffsetOut.cols; j++){
			//if (ConstOut.at<uchar>(i, j)){
			//	continue;
			//}
			OffsetOut.at<Vec2f>(i, j) = OffsetIn.at<Vec2f>(i / 2, j / 2) * 2;
			//if (HoleOut.at<uchar>(i, j)){
				int guide = GuideOut.at<uchar>(i, j);
				if (guide >= searchSpaces.size() || searchSpaces[guide].size() == 0){
					guide = 0;
				}
				Vec2f offset = OffsetOut.at<Vec2f>(i, j);
				if (i + offset[1]<0 || i + offset[1]>GuideOut.rows - 1 || j + offset[0]<0 || j + offset[0]>GuideOut.cols - 1 ||
					GuideOut.at<uchar>(i + offset[1], j + offset[0]) != guide){
					OffsetOut.at<Vec2f>(i, j) = searchSpaces[guide][rand() % searchSpaces[guide].size()];
				}
			//}
		}
	}
}
void _generateImage(Mat Offset, Mat SrcImg, Mat& DstImg, Mat Const, int patchSize, bool up){
	Mat _DstImg;
	DstImg.copyTo(_DstImg);
	float oriSize = patchSize;
	patchSize = patchSize / 2;
	Mat Count(DstImg.rows, DstImg.cols, CV_32FC1, Scalar::all(0));
	//Mat Cover(DstImg.rows, DstImg.cols, CV_32FC1, Scalar::all(0));
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
				int left = MIN(tx, patchSize), right = MIN(DstImg.cols - tx, patchSize + fix);
				int up = MIN(ty, patchSize), down = MIN(DstImg.rows - ty, patchSize + fix);
				if (sx - left<0 || sx + right>SrcImg.cols || sx - left > sx + right || sy - up<0 || sy + down>SrcImg.rows || sy - up > sy + down){
					//fprintf(stderr, "Warning:target(%d,%d) match source(%d,%d) is not valid.\n", tx, ty, sx, sy);
				}
				else{
					Mat Patch = SrcImg(Range(sy - up, sy + down), Range(sx - left, sx + right));
					Mat Dst = DstImg(Range(ty - up, ty + down), Range(tx - left, tx + right));
					Mat Cou = Count(Range(ty - up, ty + down), Range(tx - left, tx + right));
					if (Const.at<uchar>(ty, tx)) {
						Dst += Patch;
						//Cover(Range(ty - up, ty + down), Range(tx - left, tx + right)) += 1;
					}
					else {
						for (int row = 0; row < Dst.rows; ++row) {
							for (int col = 0; col < Dst.cols; ++col) {
								Dst.at<Vec3f>(row, col) += _DstImg.at<Vec3f>(ty, tx);
							}
						}
					}
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
	for (int row = 0; row < Const.rows; ++row) {
		for (int col = 0; col < Const.cols; ++col) {
			if (!Const.at<uchar>(row, col)) {
				//float ratio = 0;
				//ratio = Cover.at<float>(row, col) / (oriSize*oriSize);
				//if (ratio > 1.0) {
				//	ratio = 1.0;
				//}
				//DstImg.at<Vec3f>(row, col) = ratio * DstImg.at<Vec3f>(row, col) + (1.0 - ratio) * _DstImg.at<Vec3f>(row, col);
				DstImg.at<Vec3f>(row, col) = _DstImg.at<Vec3f>(row, col);
			}
		}
	}
}
void generateImage(Mat Offset, Mat SrcImg, Mat& DstImg, Mat Const, int patchSize, bool up){
	Mat _DstImg;
	DstImg.copyTo(_DstImg);
	float oriSize = patchSize;
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
			int coherent = -1;
			const float coherent_threshold = 1.0;
			for (int roff = -1; roff <= 1; ++roff) {
				for (int coff = -1; coff <= 1; ++coff) {
					Vec2f distance = Offset.at<Vec2f>(ty + roff, tx + coff) - offset;
					if (distance.dot(distance) <= coherent_threshold) {
						coherent++;
					}
				}
			}
			int sx = offset[0] + tx, sy = offset[1] + ty;
			int left = MIN(tx, patchSize), right = MIN(DstImg.cols - tx, patchSize + fix);
			int up = MIN(ty, patchSize), down = MIN(DstImg.rows - ty, patchSize + fix);
			if (!(sx - left<0 || sx + right>SrcImg.cols || sx - left > sx + right || sy - up<0 || sy + down>SrcImg.rows || sy - up > sy + down)) {
				Mat Patch = SrcImg(Range(sy - up, sy + down), Range(sx - left, sx + right));
				Mat Origin = _DstImg(Range(ty - up, ty + down), Range(tx - left, tx + right));
				Mat Dst = DstImg(Range(ty - up, ty + down), Range(tx - left, tx + right));
				Mat Cou = Count(Range(ty - up, ty + down), Range(tx - left, tx + right));
				Dst += Patch * float(coherent * coherent);
				Cou += coherent * coherent;
			}
		}
	}
	std::vector<Mat> channels(3);
	split(DstImg, channels);
	channels[0] /= Count;
	channels[1] /= Count;
	channels[2] /= Count;
	merge(channels, DstImg);

	//Mat Flag(Const.rows, Const.cols, CV_32FC1, Scalar::all(0));
	//for (int row = 0; row < Const.rows; ++row) {
	//	for (int col = 0; col < Const.cols; ++col) {
	//		if (Const.at<uchar>(row, col)) {
	//			for (int i = -oriSize; i < oriSize; ++i) {
	//				for (int j = -oriSize; j < oriSize; ++j) {
	//					if (row + i >= 0 && row + i < Const.rows && col + j >= 0 && col + j < Const.cols) {
	//						Flag.at<float>(row + i, col + j)++;
	//					}
	//				}
	//			}
	//		}
	//	}
	//}
	//for (int row = 0; row < Const.rows; ++row) {
	//	for (int col = 0; col < Const.cols; ++col) {
	//		if (!Const.at<uchar>(row, col)) {
	//			float ratio = (Flag.at<float>(row, col) > 0) ? 0.3 : 0;
	//			DstImg.at<Vec3f>(row, col) = ratio * DstImg.at<Vec3f>(row, col) + (1 - ratio) * _DstImg.at<Vec3f>(row, col);
	//		}
	//	}
	//}
}
void refineImage(Mat& Src, Mat& Dst, Mat& Const, int PatchSize) {
	vector<Point> boundry;
	for (int i = 0; i < Const.rows; ++i) {
		for (int j = 0; j < Const.cols; ++j) {
			if (!Const.at<uchar>(i, j)) {
				if ((i - 1 >= 0 && Const.at<uchar>(i - 1, j)) || (i + 1 < Const.rows && Const.at<uchar>(i + 1, j))
					|| (j - 1 >= 0 && Const.at<uchar>(i, j - 1)) || (j + 1 < Const.cols && Const.at<uchar>(i, j + 1))) {
					boundry.push_back(Point(i, j));
				}
			}
		}
	}
	const int RANGE = 10;
	Mat Count(Const.rows, Const.cols, CV_32FC1, Scalar::all(0));
	//Mat Canvas(Const.rows, Const.cols, CV_32FC3, Scalar::all(0));
	//for (size_t i = 0; i < boundry.size(); ++i) {
	//	int centre_row = boundry[i].x, centre_col = boundry[i].y;
	//	vector<Point> pnts;
	//	Vec3f distance(0,0,0);
	//	for (int row = MAX(0, centre_row - RANGE); row <= MIN(Const.rows - 1, centre_row + RANGE); ++row) {
	//		for (int col = MAX(0, centre_col - RANGE); col <= MIN(Const.cols - 1, centre_col + RANGE); ++col) {
	//			if (Const.at<uchar>(row, col)) {
	//				pnts.push_back(Point(row, col));
	//				distance += (Dst.at<Vec3f>(row, col) - Src.at<Vec3f>(row, col));
	//			}
	//		}
	//	}
	//	distance /= float(pnts.size());
	//	for (size_t j = 0; j < pnts.size(); ++j) {
	//		int x = pnts[j].x, y = pnts[j].y;
	//		Canvas.at<Vec3f>(x, y) += (Src.at<Vec3f>(x, y) + distance);
	//		Count.at<float>(x, y)++;
	//	}
	//}
	//for (int row = 0; row < Const.rows; ++row) {
	//	for (int col = 0; col < Const.cols; ++col) {
	//		if (Canvas.at<Vec3f>(row, col) != Vec3f(0, 0, 0)) {
	//			Dst.at<Vec3f>(row, col) = Canvas.at<Vec3f>(row, col) / Count.at<float>(row, col);
	//			for (int x = 0; x <= 2; ++x) {
	//				if (Dst.at<Vec3f>(row, col).val[x] > 1.0) {
	//					Dst.at<Vec3f>(row, col).val[x] = 1.0;
	//				}
	//			}
	//		}
	//	}
	//}
	for (size_t i = 0; i < boundry.size(); ++i) {
		int centre_row = boundry[i].x, centre_col = boundry[i].y;
		for (int row = MAX(0, centre_row - RANGE); row <= MIN(Const.rows - 1, centre_row + RANGE); ++row) {
			for (int col = MAX(0, centre_col - RANGE); col <= MIN(Const.cols - 1, centre_col + RANGE); ++col) {
				if (Const.at<uchar>(row, col)) {
					float value = 1 - float((row - centre_row) * (row - centre_row) + (col - centre_col) * (col - centre_col)) / float(2 * RANGE * RANGE);
					if (value > Count.at<float>(row, col)) {
						Count.at<float>(row, col) = value;
					}
				}
			}
		}
	}
	for (int row = 0; row < Const.rows; ++row) {
		for (int col = 0; col < Const.cols; ++col) {
			float ratio = Count.at<float>(row, col);
			Dst.at<Vec3f>(row, col) = Src.at<Vec3f>(row, col) * ratio + Dst.at<Vec3f>(row, col) * (1 - ratio);
		}
	}

	//use patchtable to refine the boundry, add source image to the image library                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
	//Mat combine = imread("combine.png");
	//combine.convertTo(combine, CV_32FC3, 1 / 255.0);
	//int rows = MAX(Src.rows, combine.rows);
	//int cols = Src.cols + combine.cols;
	//Mat library = Mat(rows, cols, CV_32FC3, Scalar::all(0));
	//Mat allow = Mat(rows, cols, CV_8UC1, Scalar::all(0));
	//for (int i = 0; i < Src.rows; ++i) {
	//	for (int j = 0; j < Src.cols; ++j) {
	//		library.at<Vec3f>(i, j) = Src.at<Vec3f>(i, j);
	//		allow.at<uchar>(i, j) = 1;
	//	}
	//}
	//for (int i = 0; i < combine.rows; ++i) {
	//	for (int j = 0; j < combine.cols; ++j) {
	//		library.at<Vec3f>(i, j + Src.cols) = combine.at<Vec3f>(i, j);
	//		allow.at<uchar>(i, j + Src.cols) = 1;
	//	}
	//}
	//PatchTable<float, float>* table = build_table(library, allow);
	//Mat Offset(Dst.rows, Dst.cols, CV_32FC2, Scalar::all(0));
	//Mat Flag(Dst.rows, Dst.cols, CV_8UC1, Scalar::all(0));
	//for (size_t i = 0; i < boundry.size(); ++i) {
	//	int centre_row = boundry[i].x, centre_col = boundry[i].y;
	//	for (int row = MAX(0, centre_row - RANGE); row <= MIN(Const.rows - 1, centre_row + RANGE); ++row) {
	//		for (int col = MAX(0, centre_col - RANGE); col <= MIN(Const.cols - 1, centre_col + RANGE); ++col) {
	//			Flag.at<uchar>(row, col) = 1;
	//		}
	//	}
	//}
	//for (int iter = 0; iter < 5; ++iter) {
	//	table_match(table, Dst, library, Offset, allow, PatchSize);
	//	generateImage(Offset, library, Dst, Flag, PatchSize);
	//}
}
void getCombination(Mat& combine, Mat& hole, Mat& guide) {
	vector<Mat> imgs;
	Mat img;
	for (int i = 1; i <= 8; ++i) {
		string name = "school";
		char num[3];
//		_itoa(i, num, 10);
        sprintf(num, "%d", i);  // Connelly: replaced Microsoft specific _itoa with sprintf
		name.append(num);
		name.append(".jpg");
		img = imread(name.c_str());
		imgs.push_back(img);
	}
	int rows = 0, cols = 0;
	for (int i = 0; i < imgs.size(); i++){
		if (imgs[i].rows > rows){
			rows = imgs[i].rows;
		}
		cols += imgs[i].cols;
	}
	hole = Mat(rows, cols, CV_8UC1);
	guide = Mat(rows, cols, CV_8UC1, Scalar::all(0));
	Mat _combine(rows, cols, CV_8UC3);

	int offset_cols = 0;
	for (int n = 0; n < imgs.size(); n++){
		int row = imgs[n].rows, col = imgs[n].cols;
		Mat img = _combine(Range(0, row), Range(offset_cols, offset_cols + col));
		for (int i = 0; i < row; i++){
			for (int j = 0; j < col; j++){
				img.at<Vec3b>(i, j) = imgs[n].at<Vec3b>(i, j);
			}
		}
		for (int i = 0; i < row - 8; i++){
			for (int j = 0; j < col - 8; j++){
				hole.at<uchar>(i, j + offset_cols) = 1;
			}
		}
		for (int i = row - 8; i < rows; i++){
			for (int j = 0; j < col; j++){
				hole.at<uchar>(i, j + offset_cols) = 0;
			}
		}
		for (int i = 0; i < row - 8; i++){
			for (int j = col - 8; j < col; j++){
				hole.at<uchar>(i, j + offset_cols) = 0;
			}
		}
		offset_cols += col;
	}
	_combine.convertTo(combine, CV_32FC3, 1 / 255.0);
	imwrite("combine.png", _combine);
}
void getReflection(vector<Mat>& imgs, vector<Mat>& holes) {
	size_t num = imgs.size();
	for (int i = 0; i < num; ++i) {
		Mat OriImg = imgs[i];
		Mat Hole = holes[i];
		Mat img1, hole1;
		img1 = Mat(OriImg.rows, OriImg.cols, OriImg.type());
		hole1 = Mat(Hole.rows, Hole.cols, Hole.type());
		for (int row = 0; row < img1.rows; ++row) {
			for (int col = 0; col < img1.cols; ++col) {
				img1.at<Vec3b>(row, col) = OriImg.at<Vec3b>(img1.rows - row - 1, col);
				hole1.at<uchar>(row, col) = Hole.at<uchar>(img1.rows - row - 1, col);
			}
		}
		imgs.push_back(img1);
		holes.push_back(hole1);
		Mat img2, hole2;
		img2 = Mat(OriImg.rows, OriImg.cols, OriImg.type());
		hole2 = Mat(Hole.rows, Hole.cols, Hole.type());
		for (int row = 0; row < img2.rows; ++row) {
			for (int col = 0; col < img2.cols; ++col) {
				img2.at<Vec3b>(row, col) = OriImg.at<Vec3b>(row, img2.cols - col - 1);
				hole2.at<uchar>(row, col) = Hole.at<uchar>(row, img2.cols - col - 1);
			}
		}
		imgs.push_back(img2);
		holes.push_back(hole2);
		Mat img3, hole3;
		img3 = Mat(OriImg.rows, OriImg.cols, OriImg.type());
		hole3 = Mat(Hole.rows, Hole.cols, Hole.type());
		for (int row = 0; row < img3.rows; ++row) {
			for (int col = 0; col < img3.cols; ++col) {
				img3.at<Vec3b>(row, col) = OriImg.at<Vec3b>(img3.rows - row - 1, img3.cols - col - 1);
				hole3.at<uchar>(row, col) = Hole.at<uchar>(img3.rows - row - 1, img3.cols - col - 1);
			}
		}
		imgs.push_back(img3);
		holes.push_back(hole3);
	}
}
IplImage* rotateImage(IplImage* img, int degree){
	double angle = degree  * CV_PI / 180.;  
	double a = sin(angle), b = cos(angle);
	int width = img->width;
	int height = img->height;
	int width_rotate = int(height * fabs(a) + width * fabs(b));
	int height_rotate = int(width * fabs(a) + height * fabs(b));
	float map[6];
	CvMat map_matrix = cvMat(2, 3, CV_32F, map);
	CvPoint2D32f center = cvPoint2D32f(width / 2, height / 2);
	cv2DRotationMatrix(center, degree, 1.0, &map_matrix);
	map[2] += (width_rotate - width) / 2;
	map[5] += (height_rotate - height) / 2;
	IplImage* img_rotate = cvCreateImage(cvSize(width_rotate, height_rotate), 8, img->nChannels);
	cvWarpAffine(img, img_rotate, &map_matrix, CV_INTER_LINEAR | CV_WARP_FILL_OUTLIERS, cvScalarAll(0));
	return img_rotate;
}
void getRotation(vector<Mat>& imgs, vector<Mat>& holes) {
	size_t num = imgs.size();
	for (int i = 0; i < num; ++i) {
		Mat OriImg = imgs[i];
		Mat Hole = holes[i];
		Mat img1 = Mat(OriImg.cols, OriImg.rows, OriImg.type());
		Mat hole1 = Mat(Hole.cols, Hole.rows, Hole.type());
		for (int row = 0; row < img1.rows; ++row) {
			for (int col = 0; col < img1.cols; ++col) {
				img1.at<Vec3b>(row, col) = OriImg.at<Vec3b>(col, row);
				hole1.at<uchar>(row, col) = Hole.at<uchar>(col, row);
			}
		}
		imgs.push_back(img1);
		holes.push_back(hole1);

		IplImage IplImg = imgs[i];
		IplImage IplHole = holes[i];
		for (int j = 1; j <= 5; ++j) {
			IplImage* _IplImg = rotateImage(&IplImg, 15 * j);
			IplImage* _IplHole = rotateImage(&IplHole, 15 * j);
			Mat _img(_IplImg, true), _hole(_IplHole, true);
			//std::cout << _img.rows << "  " << _img.cols << "  " << _hole.rows << "  " << _hole.cols << std::endl;
			//system("pause");
			imgs.push_back(_img);
			holes.push_back(_hole);
			cvReleaseImage(&_IplImg);
			cvReleaseImage(&_IplHole);
		}
	}
}
void getMultiSize(vector<Mat>& imgs, vector<Mat>& holes) {
	int num = imgs.size();
	for (size_t i = 0; i < num; ++i) {
		float rows = imgs[i].rows, cols = imgs[i].cols;
		Size size1(rows * 0.9, cols * 0.9);
		Size size2(rows * 0.8, cols * 0.8);
		Size size3(rows * 1.1, cols * 1.1);
		Size size4(rows * 1.2, cols * 1.2);
		Mat img1, img2, img3, img4;
		Mat hole1, hole2, hole3, hole4;
		resize(imgs[i], img1, size1);
		resize(imgs[i], img2, size2);
		resize(imgs[i], img3, size3);
		resize(imgs[i], img4, size4);
		resize(holes[i], hole1, size1);
		resize(holes[i], hole2, size2);
		resize(holes[i], hole3, size3);
		resize(holes[i], hole4, size4);
		imgs.push_back(img1);
		imgs.push_back(img2);
		imgs.push_back(img3);
		imgs.push_back(img4);
		holes.push_back(hole1);
		holes.push_back(hole2);
		holes.push_back(hole3);
		holes.push_back(hole4);
	}
}
void getTransformation(Mat& transformation, Mat OriImg, Mat& hole, Mat Cover, Mat& guide)
{
	vector<Mat> imgs;
	vector<Mat> holes;
	imgs.push_back(OriImg);
	holes.push_back(Cover);
	getReflection(imgs, holes);
	getRotation(imgs, holes);
	getMultiSize(imgs, holes);

	int rows = 0, cols = 0;
	for (int i = 0; i < imgs.size(); i++){
		if (imgs[i].rows > rows){
			rows = imgs[i].rows;
		}
		cols += imgs[i].cols;
	}
	hole = Mat(rows, cols, CV_8UC1);
	guide = Mat(rows, cols, CV_8UC1, Scalar::all(0));
	Mat _combine(rows, cols, CV_8UC3);

	int offset_cols = 0;
	for (int n = 0; n < imgs.size(); n++){
		int row = imgs[n].rows, col = imgs[n].cols;
		Mat img = _combine(Range(0, row), Range(offset_cols, offset_cols + col));
		for (int i = 0; i < row; i++){
			for (int j = 0; j < col; j++){
				img.at<Vec3b>(i, j) = imgs[n].at<Vec3b>(i, j);
			}
		}
		for (int i = 0; i < row - 8; i++){
			for (int j = 0; j < col - 8; j++){
				hole.at<uchar>(i, j + offset_cols) = 1 - holes[n].at<uchar>(i, j);
			}
		}
		for (int i = row - 8; i < rows; i++){
			for (int j = 0; j < col; j++){
				hole.at<uchar>(i, j + offset_cols) = 0;
			}
		}
		for (int i = 0; i < row - 8; i++){
			for (int j = col - 8; j < col; j++){
				hole.at<uchar>(i, j + offset_cols) = 0;
			}
		}
		offset_cols += col;
	}
	_combine.convertTo(transformation, CV_32FC3, 1 / 255.0);
	imwrite("hole.png", hole * 255);
	imwrite("transformation.png", _combine);
}
void assessResult(Mat& Offset, Mat& Const, float& coherence) {
	const float coherent_threshold = 1.0;
	float cnt = 0.0, sum = 0.0;
	for (int row = 0; row < Const.rows; ++row) {
		for (int col = 0; col < Const.cols; ++col) {
			if (Const.at<uchar>(row, col)) {
				for (int roff = -1; roff <= 1; ++roff) {
					for (int coff = -1; coff <= 1; ++coff) {
						if (row + roff >= 0 && row + roff < Const.rows && col + coff >= 0 && col + coff < Const.cols &&
							Const.at<uchar>(row + roff, col + coff)) {
							++sum;
							Vec2f dis = Offset.at<Vec2f>(row, col) - Offset.at<Vec2f>(row + roff, col + coff);
							if (dis.dot(dis) < coherent_threshold) {
								++cnt;
							}
						}
					}
				}
				--sum;
				--cnt;
			}
		}
	}
	coherence = cnt / sum;
}