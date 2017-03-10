#ifndef BASIC_H
#define BASIC_H

#include <cv.h>
#include <cxcore.h>
#include <highgui.h>
#include <iostream>
#include <fstream>
using namespace cv;
using namespace std;

typedef double pixel_t;//像素类型

//初始sigma
#define INIT_SIGMA 0.5
//sigma 
#define SIGMA 1.6
//intervals
#define INTERVALS 3
//r
#define RATIO 10
#define MAX_INTERPOLATION_STEPS 5 
//|D(x^)| < 0.03   0.03极值点太多
#define DXTHRESHOLD 0.04
//bins = 36
#define ORI_HIST_BINS 36    
#define ORI_SIGMA_TIMES 1.5
#define ORI_WINDOW_RADIUS 3.0 * ORI_SIGMA_TIMES 
#define ORI_SMOOTH_TIMES 2
#define ORI_PEAK_RATIO 0.8
#define FEATURE_ELEMENT_LENGTH 128
#define DESCR_HIST_BINS 8
#define DESCR_WINDOW_WIDTH 4
#define DESCR_SCALE_ADJUST 3
#define DESCR_MAG_THR 0.2
#define INT_DESCR_FCTR 512.0


#define DAt(x, y) (*(data+(y)*step+(x))) 
#define Hat(i, j) (*(H+(i)*3 + (j)))
#define At(index, x, y) (PyrAt(dog_pyr, (index), (x), (y)))
#define HIat(i, j) (*(H_inve+(i)*3 + (j)))
//抛物插值
#define Parabola_Interpolate(l, c, r) (0.5*((l)-(r))/((l)-2.0*(c)+(r))) 



struct Keypoint
{
	int octave; //关键点所在组
	int interval;// 关键点所在层

	double offset_interval;//调整后的层的增量

	int x; //x,y坐标,根据octave和interval可取的层内图像
	int y;

	//scale = sigma0*pow(2.0, o+s/S)
	double scale; //空间尺度坐标

	double dx; //特征点坐标，该坐标被缩放成原图像大小 
	double dy;

	double offset_x;
	double offset_y;

	//高斯金字塔组内各层尺度坐标，不同组的相同层的sigma值相同
	//关键点所在组的组内尺度
	double octave_scale; //offset_i;

	double ori;//方向

	int descr_length;
	double descriptor[FEATURE_ELEMENT_LENGTH]; //128 特征维数

	double val;//极值
};
#endif