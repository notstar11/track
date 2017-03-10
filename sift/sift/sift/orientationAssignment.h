#ifndef ORIENTATIONASSIGNMENT_H
#define ORIENTATIONASSIGNMENT_H

#include "basic.h"
void CalculateScale(vector<Keypoint>& features, double sigma = SIGMA, int intervals = INTERVALS);
void HalfFeatures(vector<Keypoint>& features);
bool CalcGradMagOri(const Mat& gauss, int x, int y, double& mag, double& ori);
double* CalculateOrientationHistogram(const Mat& gauss, int x, int y, int bins, int radius, double sigma);
void GaussSmoothOriHist(double *hist, int n);
double DominantDirection(double *hist, int n);
void CopyKeypoint(const Keypoint& src, Keypoint& dst);
void CalcOriFeatures(const Keypoint& keypoint, vector<Keypoint>& features, const double *hist, int n, double mag_thr);
void OrientationAssignment(vector<Keypoint>& extrema, vector<Keypoint>& features, const vector<Mat>& gauss_pyr);


#endif 