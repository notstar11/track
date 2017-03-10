#ifndef BASIC_H
#define BASIC_H

#include <cv.h>
#include <cxcore.h>
#include <highgui.h>
#include <iostream>
#include <fstream>
using namespace cv;
using namespace std;

typedef double pixel_t;//��������

//��ʼsigma
#define INIT_SIGMA 0.5
//sigma 
#define SIGMA 1.6
//intervals
#define INTERVALS 3
//r
#define RATIO 10
#define MAX_INTERPOLATION_STEPS 5 
//|D(x^)| < 0.03   0.03��ֵ��̫��
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
//�����ֵ
#define Parabola_Interpolate(l, c, r) (0.5*((l)-(r))/((l)-2.0*(c)+(r))) 



struct Keypoint
{
	int octave; //�ؼ���������
	int interval;// �ؼ������ڲ�

	double offset_interval;//������Ĳ������

	int x; //x,y����,����octave��interval��ȡ�Ĳ���ͼ��
	int y;

	//scale = sigma0*pow(2.0, o+s/S)
	double scale; //�ռ�߶�����

	double dx; //���������꣬�����걻���ų�ԭͼ���С 
	double dy;

	double offset_x;
	double offset_y;

	//��˹���������ڸ���߶����꣬��ͬ�����ͬ���sigmaֵ��ͬ
	//�ؼ�������������ڳ߶�
	double octave_scale; //offset_i;

	double ori;//����

	int descr_length;
	double descriptor[FEATURE_ELEMENT_LENGTH]; //128 ����ά��

	double val;//��ֵ
};
#endif