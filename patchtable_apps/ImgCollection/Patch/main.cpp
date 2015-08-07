#include<cstdio>
#include "PatchMatch.h"
#include "graphcuts.h"

#define index our_index     // Connelly: workaround for name conflict on Mac

//stat mat
Mat GUIImg;
Mat SrcImg;
Mat DstImg;
Mat Hole;
Mat Guide;
Mat Offset;
Mat Const;
Mat Cover;
Mat TargetImg;

//operating mat
Mat Select;

bool drawMode;
//pen;
Mat mpoint;
float r;
int pennum;
int index;

bool dragMode;
//offset
bool keep;
Vec2f start;
Vec2f offset;

//redraw the pixel(x,y) shown in the screen
void redraw(int x, int y){
	//int guide = Guide.at<uchar>(y, x);
	//if (guide > 0){
	//	switch (guide % 3){
	//	case 1:
	//		GUIImg.at<Vec3b>(y, x) = Vec3b(255, 0, 0);
	//		break;
	//	case 2:
	//		GUIImg.at<Vec3b>(y, x) = Vec3b(0, 255, 0);
	//		break;
	//	case 0:
	//		GUIImg.at<Vec3b>(y, x) = Vec3b(0, 0, 255);
	//		break;
	//	}
	//}
	//else if (Hole.at<uchar>(y, x)){
	//	GUIImg.at<Vec3b>(y, x) = Vec3b(0, 0, 0);
	//}
	//else{
	//	GUIImg.at<Vec3b>(y, x) = DstImg.at<Vec3b>(y, x);
	//}
	//if (Const.at<uchar>(y, x)){
	//	GUIImg.at<Vec3b>(y, x) = 0.5*GUIImg.at<Vec3b>(y, x) + 0.5*Vec3b(0, 255, 0);
	//}
	//if (y - offset[1] >= 0 && y - offset[1] < Select.rows && x - offset[0] >= 0 && x - offset[0] < Select.cols && Select.at<uchar>(y - offset[1], x - offset[0])){
	//	GUIImg.at<Vec3b>(y, x) = 0.5*GUIImg.at<Vec3b>(y, x) + 0.5*Vec3b(0, 0, 255);
	//}

	if (!Cover.at<uchar>(y, x)) {
		GUIImg.at<Vec3b>(y, x) = 0.5*GUIImg.at<Vec3b>(y, x) + 0.5*Vec3b(0, 255, 0);
		Const.at<uchar>(y, x) = 1;
		Cover.at<uchar>(y, x) = 1;
	}

}

void draw(int x, int y){
	for (int i = -r; i < r; i++){
		for (int j = -r; j < r; j++){
			if (i*i + j*j < r*r){
				if (x + i >= 0 && x + i < GUIImg.cols && y + j >= 0 && y + j < GUIImg.rows){
					mpoint.at<uchar>(y + j, x + i) = pennum;
					redraw(x + i, y + j);
				}
			}
		}
	}
	imshow("GUI", GUIImg);
}

void drag(){
	for (int i = 0; i < GUIImg.rows; i++){
		for (int j = 0; j < GUIImg.cols; j++){
			redraw(j, i);
		}
	}
	imshow("GUI", GUIImg);
}

void operate(){
	for (int i = 0; i < GUIImg.rows; i++){
		for (int j = 0; j < GUIImg.cols; j++){
			if (!keep && Select.at<uchar>(i, j)){
				Hole.at<uchar>(i, j) = 1;
			}
			if (i - offset[1] >= 0 && i - offset[1] < Select.rows && j - offset[0] >= 0 && j - offset[0] < Select.cols && Select.at<uchar>(i - offset[1], j - offset[0])){
				Offset.at<Vec2f>(i, j) = -offset;
				Const.at<uchar>(i, j) = 1;
				DstImg.at<Vec3b>(i, j) = SrcImg.at<Vec3b>(i - offset[1], j - offset[0]);
			}
		}
	}
	Select = Mat(SrcImg.rows, SrcImg.cols, CV_8UC1, Scalar::all(0));
	dragMode = false;
	drawMode = true;
	drag();
}

void _MouseCB(int event, int x, int y, int flags, void* param){
	static int LeftButtonPressed = false;
	if (event == CV_EVENT_LBUTTONDOWN){
		LeftButtonPressed = true;
		if (drawMode){
			draw(x, y);
		}
		if (dragMode){
			start = Vec2f(x, y) - offset;
		}
	}
	if (event == CV_EVENT_LBUTTONUP){
		LeftButtonPressed = false;
		if (dragMode){
			offset = Vec2f(x, y) - start;
			drag();
		}
	}
	if (event == CV_EVENT_MOUSEMOVE && LeftButtonPressed){
		if (drawMode){
			draw(x, y);
		}
	}
}

void _twoImagesConnection()
{
	Mat Img1 = imread("sakura4.jpg");
	Mat Img2 = imread("sakura1.jpg");
	vector<Mat> imgs;
	for (int i = 1; i <= 4; ++i) {
		string name = "sakura";
		char num[3];
		//_itoa(i, num, 10);
                sprintf(num, "%d", i);  // Connelly: replaced Microsoft specific _itoa with sprintf
		name.append(num);
		name.append(".jpg");
		Mat img = imread(name.c_str());
		imgs.push_back(img);
	}
	int optimalWidth = 0, optimalDis = INT_MAX;
	for (int w = 50; w <= 250; ++w) {
		for (size_t i = 0; i < imgs.size(); ++i) {
			for (int col = 0; col < imgs[i].cols - w; ++col) {
				int disCnt1 = 0, disCnt2 = 0;
				for (int row = 0; row < Img1.rows; ++row) {
					Vec3b dis1 = imgs[i].at<Vec3b>(row, col) - Img1.at<Vec3b>(row, Img1.cols - w);
					Vec3b dis2 = imgs[i].at<Vec3b>(row, col + w) - Img2.at<Vec3b>(row, w);
					disCnt1 += dis1.dot(dis1);
					disCnt2	+= dis2.dot(dis2);
				}
				//int disCnt = MAX(disCnt1, disCnt2);
				int disCnt = disCnt1 + disCnt2;
				if (disCnt < optimalDis) {
					optimalDis = disCnt;
					optimalWidth = w;
				}
			}
		}
	}

	const int WIDTH = optimalWidth;
	GUIImg = Mat(Img1.rows, Img1.cols + Img2.cols - WIDTH, CV_8UC3, Scalar::all(0));
	for (int row = 0; row < GUIImg.rows; ++row) {
		for (int col = 0; col < Img1.cols - WIDTH / 2; ++col) {
			GUIImg.at<Vec3b>(row, col) = Img1.at<Vec3b>(row, col);
		}
		for (int col = WIDTH / 2; col < Img2.cols; ++col) {
			GUIImg.at<Vec3b>(row, Img1.cols + col - WIDTH) = Img2.at<Vec3b>(row, col);
		}
	}
	GUIImg.copyTo(DstImg);
	GUIImg.copyTo(SrcImg);
	Const = Mat(SrcImg.rows, SrcImg.cols, CV_8UC1, Scalar::all(0));
	Cover = Mat(SrcImg.rows, SrcImg.cols, CV_8UC1, Scalar::all(0));
	Guide = Mat(SrcImg.rows, SrcImg.cols, CV_8UC1, Scalar::all(0));
	Offset = Mat(SrcImg.rows, SrcImg.cols, CV_32FC2, Scalar::all(0));
	for (int row = 0; row < GUIImg.rows; ++row) {
		for (int col = Img1.cols - WIDTH; col < Img1.cols; ++col) {
			Const.at<uchar>(row, col) = 1;
			Cover.at<uchar>(row, col) = 1;
		}
	}
	imshow("GUI", GUIImg);
}

void twoImagesConnection()
{
	int number1 = 7, number2 = 4;
	//std::cin >> number1;
	//std::cin >> number2;
	string file1("school");
	string file2("school");
	file1.append(std::to_string(number1));
	file1.append(".jpg");
	file2.append(std::to_string(number2));
	file2.append(".jpg");
	Mat Img1 = imread(file1.c_str());
	Mat Img2 = imread(file2.c_str());
	vector<Mat> imgs;
	for (int i = 1; i <= 8; ++i) {
		string name = "school";
		char num[3];
		//_itoa(i, num, 10);
                sprintf(num, "%d", i);  // Connelly: replaced Microsoft specific _itoa with sprintf
		name.append(num);
		name.append(".jpg");
		Mat img = imread(name.c_str());
		imgs.push_back(img);
	}
	int optimalWidth = 0, optimalDis = INT_MAX;
	int optimalCol = 0;
	for (int w = 50; w <= 350; ++w) {
		for (size_t i = 0; i < imgs.size(); ++i) {
			for (int col = 0; col < imgs[i].cols - w; ++col) {
				int disCnt = 0;
				for (int row = 0; row < Img1.rows; ++row) {
					Vec3b dis1 = imgs[i].at<Vec3b>(row, col) - Img1.at<Vec3b>(row, Img1.cols - 1);
					Vec3b dis2 = imgs[i].at<Vec3b>(row, col + w) - Img2.at<Vec3b>(row, 0);
					disCnt += (dis1.dot(dis1) + dis2.dot(dis2));
				}
				if (disCnt < optimalDis) {
					optimalDis = disCnt;
					optimalWidth = w;
					optimalCol = 0;
				}      
			}
		}
	}

	const int WIDTH = optimalWidth;
	GUIImg = Mat(Img1.rows, Img1.cols + Img2.cols + WIDTH, CV_8UC3, Scalar::all(0));
	for (int row = 0; row < GUIImg.rows; ++row) {
		for (int col = 0; col < Img1.cols; ++col) {
			GUIImg.at<Vec3b>(row, col) = Img1.at<Vec3b>(row, col);
		}
		for (int col = 0; col < Img2.cols; ++col) {
			GUIImg.at<Vec3b>(row, Img1.cols + col + WIDTH) = Img2.at<Vec3b>(row, col);
		}
	}
	GUIImg.copyTo(DstImg);
	GUIImg.copyTo(SrcImg);
	Const = Mat(SrcImg.rows, SrcImg.cols, CV_8UC1, Scalar::all(0));
	Cover = Mat(SrcImg.rows, SrcImg.cols, CV_8UC1, Scalar::all(0));
	Guide = Mat(SrcImg.rows, SrcImg.cols, CV_8UC1, Scalar::all(0));
	Offset = Mat(SrcImg.rows, SrcImg.cols, CV_32FC2, Scalar::all(0));
	for (int row = 0; row < GUIImg.rows; ++row) {
		for (int col = Img1.cols; col < Img1.cols + WIDTH; ++col) {
			Const.at<uchar>(row, col) = 1;
			Cover.at<uchar>(row, col) = 1;
		}
	}
	imshow("GUI", GUIImg);
}

int main(int argc, char** argv){
	SrcImg = imread("test3.jpg");
        printf("SrcImg: %dx%d\n", SrcImg.rows, SrcImg.cols);
	DstImg = SrcImg.clone();
	GUIImg = SrcImg.clone();
	Hole = Mat(SrcImg.rows, SrcImg.cols, CV_8UC1, Scalar::all(0));
	Guide = Mat(SrcImg.rows, SrcImg.cols, CV_8UC1, Scalar::all(0));
	Offset = Mat(SrcImg.rows, SrcImg.cols, CV_32FC2, Scalar::all(0));
	Const = Mat(SrcImg.rows, SrcImg.cols, CV_8UC1, Scalar::all(0));
	Cover = Mat(SrcImg.rows, SrcImg.cols, CV_8UC1, Scalar::all(0));

	Select = Mat(SrcImg.rows, SrcImg.cols, CV_8UC1, Scalar::all(0));

	drawMode = true;
	mpoint = Select;
	r = 20;
	pennum = 1;
	index = 0;

	dragMode = false;
	keep = false;
	start = Vec2f(0, 0);
	offset = Vec2f(0, 0);

        printf("Before imshow\n");
	imshow("GUI", GUIImg);
        printf("After imshow\n");
	cvSetMouseCallback("GUI", _MouseCB, NULL);
	twoImagesConnection();
	while (true){
		int key = cvWaitKey();
		switch (key){
		case '+':
			r++;
			break;
		case '-':
			if (r > 2){
				r--;
			}
			break;
		case 'g':
			mpoint = Guide;
			index++;
			pennum = index;
			break;
		case 'h':
			mpoint = Hole;
			pennum = 1;
			break;
		case 's':
			mpoint = Select;
			pennum = 1;
			break;
		case 'e':
			dragMode = false;
			drawMode = true;
			break;
		case 'f':
			mpoint = Const;
			pennum = 1;
			break;
		case 'x':
			dragMode = true;
			drawMode = false;
			keep = false;
			break;
		case 'c':
			dragMode = true;
			drawMode = false;
			keep = true;
			break;
		case 'r':
			//graphCutsBlend(DstImg, Cover);
			float maxCoherence = 0.0;
			int resIndex = 0;
			for (int idx = 0; idx < 1; ++idx) {
				float coherence = 0.0;
				DstImg.copyTo(TargetImg);
				retarget(SrcImg, TargetImg, Cover, Hole, Guide, Const, Offset, 30, 7, 100, 5, 5, 10, 10, idx, coherence);
				if (coherence > maxCoherence) {
					maxCoherence = coherence;
					resIndex = idx;
				}
			}
			imwrite("GUI.jpg", GUIImg);
			string result_name = "result_";
			result_name.append(std::to_string(resIndex));
			result_name.append("_");
			result_name.append(std::to_string(maxCoherence));
			result_name.append(".png");
			DstImg = imread(result_name.c_str());
			imwrite("result.png", DstImg);
			GUIImg = DstImg.clone();
			SrcImg = DstImg.clone();
			Hole = Mat(SrcImg.rows, SrcImg.cols, CV_8UC1, Scalar::all(0));
			Guide = Mat(SrcImg.rows, SrcImg.cols, CV_8UC1, Scalar::all(0));
			Offset = Mat(SrcImg.rows, SrcImg.cols, CV_32FC2, Scalar::all(0));
			Const = Mat(SrcImg.rows, SrcImg.cols, CV_8UC1, Scalar::all(0));
			Select = Mat(SrcImg.rows, SrcImg.cols, CV_8UC1, Scalar::all(0));
			Cover = Mat(SrcImg.rows, SrcImg.cols, CV_8UC1, Scalar::all(0));
			mpoint = Select;
			index = 0;
			dragMode = false;
			drawMode = true;
			offset = 0;
			start = 0;
			imshow("GUI", GUIImg);
			break;
		}
	}

}
