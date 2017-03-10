#ifndef DETECTEXTREME_H
#define DETECTEXTREME_H
#include "basic.h"

bool isExtremum(int x, int y, const vector<Mat>& dog_pyr, int index);
bool passEdgeResponse(int x, int y, const vector<Mat>& dog_pyr, int index, double r = RATIO);
double PyrAt(const vector<Mat>& pyr, int index, int x, int y);
void DerivativeOf3D(int x, int y, const vector<Mat>& dog_pyr, int index, double *dx);
void Hessian3D(int x, int y, const vector<Mat>& dog_pyr, int index, double *H);
bool Inverse3D(const double *H, double *H_inve);
void GetOffsetX(int x, int y, const vector<Mat>& dog_pyr, int index, double *offset_x);
double GetFabsDx(int x, int y, const vector<Mat>& dog_pyr, int index, const double* offset_x);
Keypoint* InterploationExtremum(int x, int y, const vector<Mat>& dog_pyr, int index, int octave, int interval, double dxthreshold = DXTHRESHOLD);
void DetectionLocalExtrema(const vector<Mat>& dog_pyr, vector<Keypoint>& extrema, int octaves, int intervals = INTERVALS);

#endif