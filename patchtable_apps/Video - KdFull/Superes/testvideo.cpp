
#define VIDEO_MODE 0        /* Whether to use video mode (Windows) */

#if VIDEO_MODE

#endif
#include "stdafx.h"
//#define __STDC_CONSTANT_MACROS
///*
//* Copyright (c) 2014 Stefano Sabatini
//*
//* Permission is hereby granted, free of charge, to any person obtaining a copy
//* of this software and associated documentation files (the "Software"), to deal
//* in the Software without restriction, including without limitation the rights
//* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//* copies of the Software, and to permit persons to whom the Software is
//* furnished to do so, subject to the following conditions:
//*
//* The above copyright notice and this permission notice shall be included in
//* all copies or substantial portions of the Software.
//*
//* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//* THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
//* THE SOFTWARE.
//*/

/**
* @file
* libavformat AVIOContext API example.
*
* Make libavformat demuxer access media content through a custom
* AVIOContext read callback.
* @example avio_reading.c
*/

#if VIDEO_MODE
 extern "C"
{
#include <libavcodec/avcodec.h>
#include <libavformat/avformat.h>
#include <libswscale/swscale.h>
#include <cmath>
}
#endif

#if BUILD_UNIX
#include <glob.h>
#define PATH_SEP "/"
#else
#define PATH_SEP "\\"
#endif

#include "SuperOperation.h"
#if VIDEO_MODE

#endif
#include <windows.h>
#include "time.h" 
#include <assert.h>

#if VIDEO_MODE
#include <Magick++.h>
using namespace Magick;
#endif

#define Wrong(_Expression, _Msg) (void)( (!!(_Expression)) || (_wassert(_CRT_WIDE(#_Msg), _CRT_WIDE(__FILE__), __LINE__), 0) )//Expression is what I want to get
string HighReDir = ".." PATH_SEP "high" PATH_SEP;
string HighReWholeDir = ".." PATH_SEP "high" PATH_SEP "whole" PATH_SEP;
string LowReWholeDir = ".." PATH_SEP "low" PATH_SEP;
string highResfile = "Highresolution.txt";
string LowResFile = "lowresolutionFile.txt";
string ResultFile = "result" PATH_SEP;

#if VIDEO_MODE

void SaveFrame(AVFrame *pFrame, int width, int height, int iFrame, string highReDir, ofstream& highresFile) {
	FILE *pFile;
	char szFilename[50];
	int  y;


	// Open file
	sprintf(szFilename, (highReDir + "frame%d.png").c_str(), iFrame);
	highresFile << szFilename << endl;

	pFile = fopen(szFilename, "wb");
	if (pFile == NULL)
		return;
	// Write header
	fprintf(pFile, "P6\n%d %d\n255\n", width, height);

	// Write pixel data
	for (y = 0; y < height; y++)
		fwrite(pFrame->data[0] + y*pFrame->linesize[0], 1, width * 3, pFile);

	// Close file
	fclose(pFile);
}
void RunFFMpegGetlowResVideo(char *ffmpegPath, char* inputVideo, char* outputVideo, char * size)
{
	char lpCmdLine[2000];
	sprintf(lpCmdLine, "%s -i %s - s %s %s", ffmpegPath, inputVideo, inputVideo, size, outputVideo);
	WinExec(lpCmdLine, SW_SHOW);
}
void RunFFMpegCutVideo(char *ffmpegPath, char* inputVideo, char* outputVideo, char* startime, char * lastime)
{
	char lpCmdLine[2000];
	/*char startime[50];
	char lastime[50];
	sprintf(startime, "00:00:04");
	sprintf(lastime, "00:00:03");
	char *ffmpegPath = "C:\\Users\\ll4jd\\Documents\\Visual Studio 2013\\lib set\\ffmpeg32static\\bin\\ffmpeg.exe";
	char *inputVideo = "cut.mp4";
	char *outputVideo = "bird.mp4";*/
	//sprintf(lpCmdLine, "%s -ss %s -i %s -to %s -c copy %s", ffmpegPath, st, inputVideo, kt, outputVideo);
	sprintf(lpCmdLine, "%s -ss %s -i %s -acodec copy -vcodec copy -t %s %s", ffmpegPath, startime, inputVideo, lastime, outputVideo);
	WinExec(lpCmdLine, SW_SHOW);
}
//
void CreateVideo(char *ffmpegPath)
{
	char lpCmdLine[2000];
	//sprintf(lpCmdLine, "%s -f image2 -i result\\finalresult%s.png  video.mpg", ffmpegPath, "%d");//ffmpeg -i input -c:v libx264 -preset slow -crf 22 -c:a copy output.mkv
	//ffmpeg -i finalresult%d.png -r 30 -c:v libx264 output.mp4
	sprintf(lpCmdLine, "%s -i result\\finalresult%s.png -r 30 -c:v libx264 result\\output2.mp4", ffmpegPath, "%d");
	WinExec(lpCmdLine, SW_SHOW);
}
void Crop()
{
	char lpCmdLine[2000];
	//sprintf(lpCmdLine, "%s -f image2 -i result\\finalresult%s.png  video.mpg", ffmpegPath, "%d");//ffmpeg -i input -c:v libx264 -preset slow -crf 22 -c:a copy output.mkv
	//ffmpeg -i finalresult%d.png -r 30 -c:v libx264 output.mp4
	sprintf(lpCmdLine, "mogrify -crop 200x50+128+112 low\\frame*.png");
	WinExec(lpCmdLine, SW_SHOW);
}
int ReadVideo(char* Filename, int FramePersInterval, int & HighResWidth, int &HighResCol, string highReDir)
{
	AVFormatContext *pFormatCtx = NULL;
	int             i, videoStream;
	AVCodecContext  *pCodecCtx = NULL;
	AVCodec         *pCodec = NULL;
	AVFrame         *pFrame = NULL;
	AVFrame         *pFrameRGB = NULL;
	AVPacket        packet;
	int             frameFinished;
	int             numBytes;
	uint8_t         *buffer = NULL;

	AVDictionary    *optionsDict = NULL;
	struct SwsContext      *sws_ctx = NULL;
	// Register all formats and codecs
	av_register_all();

	// Open video file
	if (avformat_open_input(&pFormatCtx, Filename, NULL, NULL) != 0)
		return -1; // Couldn't open file
	//Wrong(avformat_open_input(&pFormatCtx, Filename, NULL, NULL) == 0, "can not open file"); //Couldn't open file

	// Retrieve stream information
	if (avformat_find_stream_info(pFormatCtx, NULL) < 0)
		return -1; // Couldn't find stream information
	//Wrong(avformat_find_stream_info(pFormatCtx, NULL) >= 0, "Couldn't find stream information");// Couldn't find stream information

	// Dump information about file onto standard error
	av_dump_format(pFormatCtx, 0, Filename, 0);

	// Find the first video stream
	videoStream = -1;
	for (i = 0; i < pFormatCtx->nb_streams; i++)
		if (pFormatCtx->streams[i]->codec->codec_type == AVMEDIA_TYPE_VIDEO) {
		videoStream = i;
		break;
		}
	if (videoStream == -1)
		return -1; // Didn't find a video stream
	//Wrong(videoStream != -1, "Didn't find a video stream");// Didn't find a video stream


	// Get a pointer to the codec context for the video stream
	pCodecCtx = pFormatCtx->streams[videoStream]->codec;

	// Find the decoder for the video stream
	pCodec = avcodec_find_decoder(pCodecCtx->codec_id);
	if (pCodec == NULL) {
		fprintf(stderr, "Unsupported codec!\n");
		return -1; // Codec not found
	}
	//Wrong(pCodec != NULL, "Unsupported codec!\n");// Unsupported codec!\n

	// Open codec
	if (avcodec_open2(pCodecCtx, pCodec, &optionsDict) < 0)
		return -1; // Could not open codec
	//Wrong(avcodec_open2(pCodecCtx, pCodec, &optionsDict) >= 0, " Could not open codec!\n");//  Could not open codec

	// Allocate video frame
	pFrame = av_frame_alloc();

	// Allocate an AVFrame structure
	pFrameRGB = av_frame_alloc();
	if (pFrameRGB == NULL)
		return -1;
	//Wrong(pFrameRGB != NULL, " pFrameRGB is null!\n");//  pFrameRGB is null


	// Determine required buffer size and allocate buffer
	numBytes = avpicture_get_size(PIX_FMT_RGB24, pCodecCtx->width,
		pCodecCtx->height);
	buffer = (uint8_t *)av_malloc(numBytes*sizeof(uint8_t));

	sws_ctx =
		sws_getContext
		(
		pCodecCtx->width,
		pCodecCtx->height,
		pCodecCtx->pix_fmt,
		pCodecCtx->width,
		pCodecCtx->height,
		PIX_FMT_RGB24,
		SWS_BILINEAR,
		NULL,
		NULL,
		NULL
		);

	// Assign appropriate parts of buffer to image planes in pFrameRGB
	// Note that pFrameRGB is an AVFrame, but AVFrame is a superset
	// of AVPicture
	avpicture_fill((AVPicture *)pFrameRGB, buffer, PIX_FMT_RGB24,
		pCodecCtx->width, pCodecCtx->height);

	// Read frames and save first five frames to disk
	i = 0;
	int TotalFrameNumber = 0;
	ofstream highresFile(highResfile);

	while (av_read_frame(pFormatCtx, &packet) >= 0) {
		// Is this a packet from the video stream?
		if (packet.stream_index == videoStream) {
			// Decode video frame
			avcodec_decode_video2(pCodecCtx, pFrame, &frameFinished,
				&packet);

			// Did we get a video frame?
			if (frameFinished) {
				// Convert the image from its native format to RGB
				sws_scale
					(
					sws_ctx,
					(uint8_t const * const *)pFrame->data,
					pFrame->linesize,
					0,
					pCodecCtx->height,
					pFrameRGB->data,
					pFrameRGB->linesize
					);

				// Save the frame to disk
				if ((i) % FramePersInterval == 0 || i == 0)
				{
					SaveFrame(pFrameRGB, pCodecCtx->width, pCodecCtx->height, i, highReDir, highresFile);
					HighResWidth = pCodecCtx->width;
					HighResCol = pCodecCtx->height;
				}
				i++;
			}
		}

		// Free the packet that was allocated by av_read_frame
		av_free_packet(&packet);
	}
	highresFile.close();
	cout << "TotalNumber is:" << i << endl;
	// Free the RGB image
	av_free(buffer);
	av_free(pFrameRGB);

	// Free the YUV frame
	av_free(pFrame);

	// Close the codec
	avcodec_close(pCodecCtx);

	// Close the video file
	avformat_close_input(&pFormatCtx);
	return i;
}
#endif

#define FramePersInterval 1
int DownSampleFactor = 1;
int FactorHandL = 8;
int DownSampleFactorForLow = DownSampleFactor*FactorHandL;
int getAllFilesNumber(string path, vector<string>& files)
{
#if BUILD_UNIX
    int fileNumber = 0;
    glob_t glob_result;
    string pattern = path + "/*.png";
    glob(pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
    for (int i = 0; i < glob_result.gl_pathc; i++) {
        files.push_back(string(glob_result.gl_pathv[i]));
        fileNumber++;
    }
    printf("fileNumber=%d\n", fileNumber);
    return fileNumber;
#else
    //文件句柄
	long   hFile = 0;
	int fileNumber = 0;
	//文件信息  
	struct _finddata_t fileinfo;  //很少用的文件信息读取结构
	string p;  //string类很有意思的一个赋值函数:assign()，有很多重载版本
	if ((hFile = _findfirst(p.assign(path).append("\\*").c_str(), &fileinfo)) != -1)
	{
		do
		{
			if ((fileinfo.attrib &  _A_SUBDIR))  //比较文件类型是否是文件夹
			{
				if (strcmp(fileinfo.name, ".") != 0 && strcmp(fileinfo.name, "..") != 0)
				{
					//files.push_back(p.assign(path).append("\\").append(fileinfo.name));
					//getAllFiles(p.assign(path).append("\\").append(fileinfo.name), files);
				}
			}
			else
			{
				files.push_back(p.assign(path).append("\\").append(fileinfo.name));
				fileNumber++;
			}
		} while (_findnext(hFile, &fileinfo) == 0);  //寻找下一个，成功返回0，否则-1
		_findclose(hFile);
	}
	return fileNumber;
#endif
}
float ComputeDifference(cv::Mat OrigHigh, cv::Mat Test)
{
	float diff = 0;
	for (int r = 0; r < OrigHigh.rows; r++)
	{
		for (int c = 0; c < OrigHigh.cols; c++)
		{
			diff += (OrigHigh.at<float>(r, c) - Test.at<float>(r, c))*(OrigHigh.at<float>(r, c) - Test.at<float>(r, c));
		}
	}
	diff = diff / (OrigHigh.rows*OrigHigh.cols);
	return diff;
}
float ComputeDifference_psnr(cv::Mat OrigHigh, cv::Mat Test)
{
	float diff = 0;
	for (int n = 0; n < 3; n++)
	{
		for (int r = 0; r < OrigHigh.rows; r++)
		{
			for (int c = 0; c < OrigHigh.cols; c++)
			{
				diff += (OrigHigh.at<cv::Vec3f>(r, c)[n] - Test.at<cv::Vec3f>(r, c)[n])*(OrigHigh.at<cv::Vec3f>(r, c)[n] - Test.at<cv::Vec3f>(r, c)[n]);
			}
		}
	}
	diff = diff / (OrigHigh.rows*OrigHigh.cols*3);
	float a = ((255.0)*(255.0)) / diff;
	float b = log10(a) * 10;
	return b;
}
void ComputeDiff_Patchtable()
{
	string OrighHigh = "result" PATH_SEP "high";

	vector<string> Highfiles;
	vector<string> Patchfiles;
	int filen_Orign = getAllFilesNumber(OrighHigh.c_str(), Highfiles);
	for (int n = 0; n < 25; n++)
	{
		string PatchTable = "result" PATH_SEP "pca20_ndims10_kco40_prop_iters4_i3";
		char c[20];
		sprintf(c, "finalresult%d", n);
		Patchfiles.push_back(PatchTable.append(PATH_SEP).append(c).append(".png"));
	}
	vector< cv::Mat> LumYcYb;
	ofstream LowResfile("error.txt", std::ofstream::out | std::ofstream::app);
	for (int n = 0; n < filen_Orign; n++)
	{
		/* cv::Mat im_HighImg = cv::imread(Highfiles[n], 1);
		 im_HighImg.convertTo(im_HighImg, CV_32FC3);
		 cv::cvtColor(im_HighImg, im_HighImg, cv::COLOR_BGR2YCrCb);
		 split(im_HighImg, LumYcYb);
		 cv::Mat im_HighLum = LumYcYb[0].clone();
		 cv::Mat im_HighCR = LumYcYb[1].clone();
		 cv::Mat im_HighCB = LumYcYb[2].clone();
		 LumYcYb.clear();*/
		cv::Mat im_HighImg = cv::imread(Highfiles[n], CV_LOAD_IMAGE_COLOR);
		im_HighImg.convertTo(im_HighImg, CV_32FC3);
		//cv::cvtColor(im_HighImg, im_HighImg, cv::COLOR_BGR2Lab);

		/*cv::Mat im_PatchImg = cv::imread(Patchfiles[n], 1);
		im_PatchImg.convertTo(im_PatchImg, CV_32FC3);
		cv::cvtColor(im_PatchImg, im_PatchImg, cv::COLOR_BGR2YCrCb);
		split(im_PatchImg, LumYcYb);
		cv::Mat im_PatchLum = LumYcYb[0].clone();
		cv::Mat im_PatchCR = LumYcYb[1].clone();
		cv::Mat im_PatchCB = LumYcYb[2].clone();
		LumYcYb.clear();*/
		//float error = ComputeDifference(im_HighLum, im_PatchLum);
		cv::Mat im_PatchImg = cv::imread(Patchfiles[n], CV_LOAD_IMAGE_COLOR);
		im_PatchImg.convertTo(im_PatchImg, CV_32FC3);
		//cv::cvtColor(im_PatchImg, im_PatchImg, cv::COLOR_BGR2Lab);
		float error = ComputeDifference_psnr(im_HighImg, im_PatchImg);
		LowResfile << error << endl;
	}
	LowResfile.close();
}
void ComputeDiff_highkdtree()
{
	string OrighHigh = ".." PATH_SEP "high";

	vector<string> Highfiles;
	vector<string> Patchfiles;
	int filen_Orign = getAllFilesNumber(OrighHigh.c_str(), Highfiles);

	for (int n = 0; n < filen_Orign; n++)
	{
		string PatchTable = "result";
		char c[20];
		sprintf(c, "finalresult%d", n);
		Patchfiles.push_back(PatchTable.append(PATH_SEP).append(c).append(".png"));
	}
	vector< cv::Mat> LumYcYb;
	ofstream LowResfile("error.txt", std::ofstream::out | std::ofstream::app);
	for (int n = 0; n < filen_Orign; n++)
	{
		/*cv::Mat im_HighImg = cv::imread(Highfiles[n], 1);
		im_HighImg.convertTo(im_HighImg, CV_32FC3);
		cv::cvtColor(im_HighImg, im_HighImg, cv::COLOR_BGR2YCrCb);
		split(im_HighImg, LumYcYb);
		cv::Mat im_HighLum = LumYcYb[0].clone();
		cv::Mat im_HighCR = LumYcYb[1].clone();
		cv::Mat im_HighCB = LumYcYb[2].clone();
		LumYcYb.clear();*/
		cv::Mat im_HighImg = cv::imread(Highfiles[n], CV_LOAD_IMAGE_COLOR);
		im_HighImg.convertTo(im_HighImg, CV_32FC3);

		/*cv::Mat im_PatchImg = cv::imread(Patchfiles[n], 1);
		im_PatchImg.convertTo(im_PatchImg, CV_32FC3);
		cv::cvtColor(im_PatchImg, im_PatchImg, cv::COLOR_BGR2YCrCb);
		split(im_PatchImg, LumYcYb);
		cv::Mat im_PatchLum = LumYcYb[0].clone();
		cv::Mat im_PatchCR = LumYcYb[1].clone();
		cv::Mat im_PatchCB = LumYcYb[2].clone();
		LumYcYb.clear();*/
		cv::Mat im_PatchImg = cv::imread(Patchfiles[n], CV_LOAD_IMAGE_COLOR);
		im_PatchImg.convertTo(im_PatchImg, CV_32FC3);
		/*float error = ComputeDifference(im_HighLum, im_PatchLum);*/
		float error = ComputeDifference_psnr(im_HighImg, im_PatchImg);
		LowResfile << error << endl;
	}
	LowResfile.close();
}
void wirte()
{
	cv::Mat test(10, 10, CV_32FC1);
	cv::randu(test, 0, 1);
	std::cout << test << std::endl;
	assert(test.isContinuous());
	//TrainFeatureData = cv::Mat::zeros(KeyRows, im_PCADIMS, CV_32FC1);
	std::ofstream myFile1("data1.bin", std::ios::out | std::ios::binary);
	for (int r = 0; r < 10; r++)
	{
		myFile1.write(reinterpret_cast<char*>(test.ptr(r)), test.cols*test.elemSize());
	}
		myFile1.close();
	float buf[100];
	std::ifstream myFile("data1.bin", std::ios::in | std::ios::binary);
	myFile.read(reinterpret_cast<char*>(buf), sizeof(buf));
    myFile.close();
    cv::Mat test1(10, 10, CV_32FC1, buf);
	std::cout << test1 << std::endl;
	cv::waitKey();
}
int main(int argc, char *argv[]) {

	char startime[50];
	char lastime[50];
	char size[50];
	sprintf(startime, "00:00:03");//start time 3D rendered：04
	sprintf(lastime, "00:00:09");//how long  3D rendered：03

	char *ffmpegPath = "C:\\Users\\ll4jd\\Documents\\Visual Studio 2013\\lib set\\ffmpeg32static\\bin\\ffmpeg.exe";
	char *inputVideo = "original.mp4";
	char *HighResVideo = "birdhigh.mp4";
	int HighResWidth = 0, HighResCol = 0;

	//RunFFMpegCutVideo(ffmpegPath, inputVideo, HighResVideo, startime, lastime);//step1: Cut the approriate video named outputVideo;
	//int lowNumber = ReadVideo(HighResVideo, 1, HighResWidth, HighResCol, HighReWholeDir);//step2: Store the highresolution pictures; how many pictures per second//all the highresolution pics prepared for the lowresolution
	//int highNumber = ReadVideo(HighResVideo, FramePersInterval, HighResWidth, HighResCol, HighReDir);//the highresolution samples
	//int LowResWidth = HighResWidth / DownSampleFactorForLow;
	//int LowResCols = HighResCol / DownSampleFactorForLow;
	////cout << LowResWidth << " " << LowResCols << " " << highNumber << " " << lowNumber;
	//sprintf(size, "%dx%d", LowResWidth, LowResCols);


	//InitializeMagick(*argv);//step3: Get lowresolution part. //RUN this code lonely and active the number 
	//ofstream LowResfile(LowResFile);
	////int lowNumber = 184;
	////just run the following code when used
	///*int highNumber = 184;
	//int LowResWidth = 444;
	//int LowResCols = 250;*/
	//char name[50];
	//Image master;
	//for (int n = 0; n < lowNumber; n++)
	//{
	//	sprintf(name, (HighReWholeDir + "frame%d.png").c_str(), n);
	//	master.read(name);
	//	Magick::FilterTypes filter(LanczosFilter);
	//	master.filterType(filter);
	//	Geometry geometry(LowResWidth, LowResCols);
	//	master.zoom(geometry);
	//	sprintf(name, (LowReWholeDir + "frame%d.png").c_str(), n);
	//	master.write(name);
	//	LowResfile << name << endl;
	//	//step4: crop small highresolution image   convert - crop 250x400 + 850 + 550 frame0.png fram0.png or mogrify -crop 672x512+2240+800 frame*.png
	//	//mogrify -crop 200x50+128+112 frame*.png
	//	/*sprintf(name, (LowReWholeDir + "frame%d.png").c_str(), n);
	//	master.read(name);
	//	Geometry geometry2(128, 112, 200, 50);
	//	master.crop(geometry2);
	//	sprintf(name, (LowReWholeDir + "frame%d.png").c_str(), n);
	//	master.write(name);
	//	LowResfile << name<<endl;*/
	//}
	//cout << LowResWidth << " " << LowResCols;
	//LowResfile.close();

	///************************/
	//////this is part is used for change the size of training pics 
	//////this step must be added because if the training set too big then it will be crushed
	///************************/
	////int lowNumber = 200;//Just test this only part
	////int highNumber = 200;
	////int LowResWidth = 256;
	////int LowResCols = 144;
	////HighResWidth = 1280;
	////HighResCol = 720;
	////InitializeMagick(*argv);
	////char name[50];
	////Image master;

	//for (int n = 0; n < highNumber / FramePersInterval; n++)
	//{
	//	sprintf(name, (HighReDir + "frame%d.png").c_str(), FramePersInterval*n);
	//	master.read(name);
	//	Magick::FilterTypes filter(LanczosFilter);
	//	master.filterType(filter);
	//	Geometry geometry(HighResWidth / DownSampleFactor, HighResCol / DownSampleFactor);
	//	//cout << HighResCol << " " << DownSampleFactor;
	//	master.zoom(geometry);
	//	sprintf(name, (HighReDir + "frame%d.png").c_str(), FramePersInterval*n);
	//	master.write(name);
	//}

	SuperOperation *s = new SuperOperation(highResfile, LowResFile, 9, "Highresultion.png", FramePersInterval);//step5: Our algorithm begin. 
	s->ProcessImage(ResultFile);
	delete s;

	//wirte();
	//ComputeDiff_Patchtable();
	ComputeDiff_highkdtree();
	//ComputeDiff_Patchtable();
	//char name1[15]; char name2[15];
	//step7: Create video
	//method 1:
	//CreateVideo(ffmpegPath);//This function can work but because the png is RGB channel and didnot convert to yuv channel. It is just so-so
	//method 2: to create myData.h264 It's very poor!
	//Step8: cubic pic:
	//string OrighHigh = "low";
	//char name[50];
	//vector<string> Highfiles;
	//int filen_Orign = getAllFilesNumber(OrighHigh.c_str(), Highfiles);
	//for (int n = 0; n < filen_Orign; n++)
	//{
	//	cv::Mat original=cv::imread(Highfiles[n]);
	//	cv::Mat out;
	//	cv::resize(original, out, cv::Size(4284, 2844), 0, 0, cv::INTER_CUBIC);
	//	//cout << HighResCol << " " << DownSampleFactor;
	//	
	//	sprintf(name, (OrighHigh + "frame%d.png").c_str(), n);
	//	cv::imwrite(name, out);
	//}
}
