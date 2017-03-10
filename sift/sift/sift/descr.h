#ifndef DESCR_H
#define DESCR_H
#include "basic.h"
#include "orientationAssignment.h"


void InterpHistEntry(double ***hist, double xbin, double ybin, double obin, double mag, int bins, int d);
double*** CalculateDescrHist(const Mat& gauss, int x, int y, double octave_scale, double ori, int bins, int width);
void NormalizeDescr(Keypoint& feat);
void HistToDescriptor(double ***hist, int width, int bins, Keypoint& feature);
void DescriptorRepresentation(vector<Keypoint>& features, const vector<Mat>& gauss_pyr, int bins, int width);

#endif