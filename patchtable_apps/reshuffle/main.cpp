#include<cstdio>
#include"reshuffle.h"
#define index main_index	/* Workaround for linker conflict on Mac */

//state mat
Mat GUIImg;
Mat SrcImg;
Mat DstImg;
Mat Hole;
Mat Guide;
Mat Offset;
Mat Const;

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
	int guide = Guide.at<uchar>(y, x);
	if (guide > 0){
		switch (guide % 3){
		case 1:
			GUIImg.at<Vec3b>(y, x) = Vec3b(255, 0, 0);
			break;
		case 2:
			GUIImg.at<Vec3b>(y, x) = Vec3b(0, 255, 0);
			break;
		case 0:
			GUIImg.at<Vec3b>(y, x) = Vec3b(0, 0, 255);
			break;
		}
	}
	else if (Hole.at<uchar>(y, x)){
		GUIImg.at<Vec3b>(y, x) = Vec3b(0, 0, 0);
	}
	else{
		GUIImg.at<Vec3b>(y, x) = DstImg.at<Vec3b>(y, x);
	}
	if (Const.at<uchar>(y, x)){
		GUIImg.at<Vec3b>(y, x) = 0.5*GUIImg.at<Vec3b>(y, x) + 0.5*Vec3b(0, 255, 0);
	}
	if (y - offset[1] >= 0 && y - offset[1] < Select.rows && x - offset[0] >= 0 && x - offset[0] < Select.cols && Select.at<uchar>(y - offset[1], x - offset[0])){
		GUIImg.at<Vec3b>(y, x) = 0.5*GUIImg.at<Vec3b>(y, x) + 0.5*Vec3b(0, 0, 255);
	}
}

void draw(int x, int y){
	for (int i = -r; i < r; i++){
		for (int j = -r; j < r; j++){
			if (i*i + j*j < r*r){
				if (x + i >= 0 && x + i < mpoint.cols && y + j >= 0 && y + j < mpoint.rows){
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

int main(int argc, char** argv){
	SrcImg = imread("Image.png");
	DstImg = SrcImg.clone();
	GUIImg = SrcImg.clone();
	Hole = Mat(SrcImg.rows, SrcImg.cols, CV_8UC1, Scalar::all(0));
	Guide = Mat(SrcImg.rows, SrcImg.cols, CV_8UC1, Scalar::all(0));
	Offset = Mat(SrcImg.rows, SrcImg.cols, CV_32FC2, Scalar::all(0));
	Const = Mat(SrcImg.rows, SrcImg.cols, CV_8UC1, Scalar::all(0));

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

	imshow("GUI", GUIImg);
	cvSetMouseCallback("GUI", _MouseCB, NULL);
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
			if (offset != Vec2f(0, 0)){
				operate();
			}
			retarget(SrcImg, DstImg, Hole, Guide, Const, Offset, 100, 7, 100, 5, 5);
			GUIImg = DstImg.clone();
			SrcImg = DstImg.clone();
			Hole = Mat(SrcImg.rows, SrcImg.cols, CV_8UC1, Scalar::all(0));
			Guide = Mat(SrcImg.rows, SrcImg.cols, CV_8UC1, Scalar::all(0));
			Offset = Mat(SrcImg.rows, SrcImg.cols, CV_32FC2, Scalar::all(0));
			Const = Mat(SrcImg.rows, SrcImg.cols, CV_8UC1, Scalar::all(0));
			Select = Mat(SrcImg.rows, SrcImg.cols, CV_8UC1, Scalar::all(0));
			mpoint = Select;
			index = 0;
			dragMode = false;
			drawMode = true;
			offset = 0;
			start = 0;
			imshow("GUI", GUIImg);
		}

	}

}
