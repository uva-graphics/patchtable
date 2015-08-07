#include "graphcuts.h"
#include <map>

void match_init(Mat& Dst0, Mat& Src0, Mat& Cover0, int& rowOff, int& colOff)
{
	Mat Dst, Src, Cover;
	pyrDown(Dst0, Dst);
	pyrDown(Src0, Src);
	pyrDown(Cover0, Cover);
	for (int i = 0; i < 2; ++i) {
		pyrDown(Dst, Dst);
		pyrDown(Src, Src);
		pyrDown(Cover, Cover);
	}

	int minRow = Cover.rows, minCol = Cover.cols, maxRow = 0, maxCol = 0;
	vector<Point> boundry;
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
		}
	}

	Point loc;
	float minDis = 1e6;
	for (int row = 0; row < Src.rows - maxRow + minRow; ++row) {
		for (int col = 0; col < Src.cols - maxCol + minCol; ++col) {
			int rowOff = row - minRow, colOff = col - minCol;
			float dis = 0.0;
			for (size_t i = 0; i < boundry.size(); ++i) {
				int x = boundry[i].x, y = boundry[i].y;
				Vec3f diff = Dst.at<Vec3f>(x, y) - Src.at<Vec3f>(x + rowOff, y + colOff);
				dis += diff.dot(diff);
			}
			if (dis < minDis) {
				minDis = dis;
				loc.x = row; loc.y = col;
			}
		}
	}
	rowOff = (loc.x - minRow) * 8;
	colOff = (loc.y - minCol) * 8;
	for (int i = 0; i < Dst0.rows; ++i) {
		for (int j = 0; j < Dst0.cols; ++j) {
			if (Cover0.at<uchar>(i, j)) {
				Dst0.at<Vec3f>(i, j) = Src0.at<Vec3f>(i + rowOff, j + colOff);
			}
		}
	}
}

void getCombinationImage(Mat& combine, Mat& hole) {
	vector<Mat> imgs;
	Mat img;
	for (int i = 1; i <= 20; ++i) {
		string name = "f";
		char num[3];
		//_itoa(i, num, 10);
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
	hole = Mat(rows, cols, CV_8UC1, Scalar::all(0));
	Mat _combine(rows, cols, CV_8UC3, Scalar::all(0));
	combine = Mat(rows, cols, CV_32FC3, Scalar::all(0));

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

	_combine.convertTo(combine, CV_32FC3, 1.0 / 255.0);
	imwrite("combine.png", _combine);
}

void optimal_boundry(Mat& Dst, Mat& Ori, Mat& Cover, Mat& Select)
{
	Cover.copyTo(Select);
	Mat Flag(Cover.rows, Cover.cols, CV_8UC1, Scalar::all(0));
	vector<Point> outer_boundry;
	for (int i = 0; i < Cover.rows; ++i) {
		for (int j = 0; j < Cover.cols; ++j) {
			if (Cover.at<uchar>(i, j)) {
				if ((i - 1 >= 0 && !Cover.at<uchar>(i - 1, j)) || (i + 1 < Cover.rows && !Cover.at<uchar>(i + 1, j))
					|| (j - 1 >= 0 && !Cover.at<uchar>(i, j - 1)) || (j + 1 < Cover.cols && !Cover.at<uchar>(i, j + 1))) {
					outer_boundry.push_back(Point(i, j));
				}
			}
		}
	}
	const int RANGE = 10;
	std::map<int, int> pnts;
	int num = 0;
	for (size_t i = 0; i < outer_boundry.size(); ++i) {
		int bx = outer_boundry[i].x, by = outer_boundry[i].y;
		for (int x = MAX(0, bx - RANGE); x <= MIN(Flag.rows - 1, bx + RANGE); ++x) {
			for (int y = MAX(0, by - RANGE); y <= MIN(Flag.cols - 1, by + RANGE); ++y) {
				if (Cover.at<uchar>(x, y)) {
					int p = x * Flag.cols + y;
					if (pnts.find(p) == pnts.end()) {
						Flag.at<uchar>(x, y) = 1;
						pnts[p] = num;
						++num;
					}
				}
			}
		}
	}
	vector<Point> inner_boundry;
	for (int i = 0; i < Flag.rows; ++i) {
		for (int j = 0; j < Flag.cols; ++j) {
			if (Flag.at<uchar>(i, j)) {
				if ((i - 1 >= 0 && Cover.at<uchar>(i - 1, j) && !Flag.at<uchar>(i - 1, j)) || (i + 1 < Cover.rows && Cover.at<uchar>(i + 1, j) && !Flag.at<uchar>(i + 1, j))
					|| (j - 1 >= 0 && Cover.at<uchar>(i, j - 1) && !Flag.at<uchar>(i, j - 1)) || (j + 1 < Cover.cols && Cover.at<uchar>(i, j + 1) && !Flag.at<uchar>(i, j + 1))) {
					inner_boundry.push_back(Point(i, j));
				}
			}
		}
	}

	Graph<float, float, float> graph(num, num * 4);
	graph.add_node(num);
	for (size_t i = 0; i < outer_boundry.size(); ++i) {
		graph.add_tweights(pnts[outer_boundry[i].x * Flag.cols + outer_boundry[i].y], 10000, 0);
	}
	for (size_t i = 0; i < inner_boundry.size(); ++i) {
		graph.add_tweights(pnts[inner_boundry[i].x * Flag.cols + inner_boundry[i].y], 0, 10000);
	}
	Mat Edge(num, num, CV_8UC1, Scalar::all(0));
	std::map<int, int>::iterator it = pnts.begin();
	while (it != pnts.end()) {
		int x = it->first / Flag.cols, y = it->first % Flag.cols;
		int idx1 = it->second;
		if (x > 0 && Flag.at<uchar>(x - 1, y)) {
			int idx2 = pnts[(x - 1) * Flag.cols + y];
			if (!Edge.at<uchar>(idx1, idx2)) {
				Edge.at<uchar>(idx1, idx2) = 1;
				Edge.at<uchar>(idx2, idx1) = 1;
				Vec3f dis1 = Dst.at<Vec3f>(x, y) - Ori.at<Vec3f>(x, y);
				Vec3f dis2 = Dst.at<Vec3f>(x - 1, y) - Ori.at<Vec3f>(x - 1, y);
				float cap = sqrt(dis1.dot(dis1)) + sqrt(dis2.dot(dis2));
				graph.add_edge(idx1, idx2, cap, cap);
			}
		}
		if (y > 0 && Flag.at<uchar>(x, y - 1)) {
			int idx2 = pnts[x * Flag.cols + y - 1];
			if (!Edge.at<uchar>(idx1, idx2)) {
				Edge.at<uchar>(idx1, idx2) = 1;
				Edge.at<uchar>(idx2, idx1) = 1;
				Vec3f dis1 = Dst.at<Vec3f>(x, y) - Ori.at<Vec3f>(x, y);
				Vec3f dis2 = Dst.at<Vec3f>(x, y - 1) - Ori.at<Vec3f>(x, y - 1);
				float cap = sqrt(dis1.dot(dis1)) + sqrt(dis2.dot(dis2));
				graph.add_edge(idx1, idx2, cap, cap);
			}
		}
		if (x < Flag.rows - 1 && Flag.at<uchar>(x + 1, y)) {
			int idx2 = pnts[(x + 1) * Flag.cols + y];
			if (!Edge.at<uchar>(idx1, idx2)) {
				Edge.at<uchar>(idx1, idx2) = 1;
				Edge.at<uchar>(idx2, idx1) = 1;
				Vec3f dis1 = Dst.at<Vec3f>(x, y) - Ori.at<Vec3f>(x, y);
				Vec3f dis2 = Dst.at<Vec3f>(x + 1, y) - Ori.at<Vec3f>(x + 1, y);
				float cap = sqrt(dis1.dot(dis1)) + sqrt(dis2.dot(dis2));
				graph.add_edge(idx1, idx2, cap, cap);
			}
		}
		if (y < Flag.cols - 1 && Flag.at<uchar>(x, y + 1)) {
			int idx2 = pnts[x * Flag.cols + y + 1];
			if (!Edge.at<uchar>(idx1, idx2)) {
				Edge.at<uchar>(idx1, idx2) = 1;
				Edge.at<uchar>(idx2, idx1) = 1;
				Vec3f dis1 = Dst.at<Vec3f>(x, y) - Ori.at<Vec3f>(x, y);
				Vec3f dis2 = Dst.at<Vec3f>(x, y + 1) - Ori.at<Vec3f>(x, y + 1);
				float cap = sqrt(dis1.dot(dis1)) + sqrt(dis2.dot(dis2));
				graph.add_edge(idx1, idx2, cap, cap);
			}
		}
		++it;
	}
	graph.maxflow();
	it = pnts.begin();
	while (it != pnts.end()) {
		int x = it->first / Flag.cols, y = it->first % Flag.cols;
		int idx = it->second;
		if (graph.what_segment(idx) == graph.SOURCE) {
			Select.at<uchar>(x, y) = 0;
			Dst.at<Vec3f>(x, y) = Ori.at<Vec3f>(x, y);
		}
		++it;
	}
}

bool is_boundry(Mat& Cover, int row, int col)
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

void poisson_blending(Mat& Dst, Mat& Ori, Mat& Cover, Mat& Select)
{
	vector<Point> boundry;
	bool start = false;
	for (int row = 0; row < Select.rows; ++row) {
		for (int col = 0; col < Select.cols; ++col) {
			if (is_boundry(Select, row, col)) {
				boundry.push_back(Point(row, col));
				start = true;
				break;
			}
		}
		if (start) {
			break;
		}
	}
	std::cout << "start to get the boundry." << std::endl;
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
		for (int x = MIN(row + 1, Select.rows - 1); x >= MAX(row - 1, 0); --x) {
			for (int y = MIN(col + 1, Select.cols - 1); y >= MAX(col - 1, 0); --y) {
				if ((x != row || y != col) && (x != _row || y != _col)) {
					if (is_boundry(Select, x, y)) {
						boundry.push_back(Point(x, y));
						std::cout << x << "  " << y << std::endl;
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
	std::cout << "find boundry." << std::endl;
	vector<Point> cover_pnts;
	for (int i = 0; i < Select.rows; ++i) {
		for (int j = 0; j < Select.cols; ++j) {
			if (Select.at<uchar>(i, j)) {
				cover_pnts.push_back(Point(i, j));
			}
		}
	}
	vector<Vec3f> boundry_offset;
	for (size_t i = 0; i < boundry.size(); ++i) {
		int x = boundry[i].x, y = boundry[i].y;
		Vec3f offset = Ori.at<Vec3f>(x, y) - Dst.at<Vec3f>(x, y);
		boundry_offset.push_back(offset);
	}

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
			color += (weights[j] / weight_sum) * boundry_offset[i];
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
	for (int i = 0; i < Cover.rows; ++i) {
		for (int j = 0; j < Cover.cols; ++j) {
			if (Cover.at<uchar>(i, j) && !Select.at<uchar>(i, j)) {
				Dst.at<Vec3f>(i, j) = Ori.at<Vec3f>(i, j);
			}
		}
	}
}

void graphCutsBlend(Mat& DstImg, Mat& Cover)
{
	Mat Src, Hole;
	Mat _DstImg(DstImg.rows, DstImg.cols, CV_32FC3);
	Mat OriImg(DstImg.rows, DstImg.cols, CV_32FC3);
	DstImg.convertTo(_DstImg, CV_32FC3, 1.0 / 255.0);
	_DstImg.copyTo(OriImg);
	int rowOff, colOff;
	Mat Select(DstImg.rows, DstImg.cols, CV_8UC1);

	getCombinationImage(Src, Hole);
	match_init(_DstImg, Src, Cover, rowOff, colOff);
	optimal_boundry(_DstImg, OriImg, Cover, Select);
	std::cout << "graph cuts over." << std::endl; 
	//poisson_blending(_DstImg, OriImg, Cover, Select);

	Mat resultImg(DstImg.rows, DstImg.cols, CV_8UC3);
	_DstImg.convertTo(resultImg, CV_8UC3, 255);
	imwrite("result_graphcuts.png", resultImg);
}
