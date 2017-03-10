#ifndef CREATEPYRAMID_H
#define CREATEPYRAMID_H

#include "basic.h"
#include "preprocess.h"

void CreateInitSmoothGray(const Mat &src, Mat &dst, double sigma = SIGMA);
void GaussianPyramid(const Mat &src, vector<Mat>&gauss_pyr, int octaves, int intervals = INTERVALS, double sigma = SIGMA);
void Sub(const Mat& a, const Mat& b, Mat & c);
void DogPyramid(const vector<Mat>& gauss_pyr, vector<Mat>& dog_pyr, int octaves, int intervals = INTERVALS);

#endif