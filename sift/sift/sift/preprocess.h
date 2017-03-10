#ifndef PREPROCESS_H
#define PREPROCESS_H
#include "basic.h"

void ConvertToGray(const Mat& src, Mat& dst);
void DownSample(const Mat& src, Mat& dst);
void UpSample(const Mat &src, Mat &dst);
void GaussianTemplateSmooth(const Mat &src, Mat &dst, double sigma);
void GaussianSmooth2D(const Mat &src, Mat &dst, double sigma);
void GaussianSmooth(const Mat &src, Mat &dst, double sigma);

#endif