#ifndef SIFT_H
#define SIFT_H

#include "createPyramid.h"
#include "descr.h"
#include "detectExtreme.h"
#include "preprocess.h"

bool FeatureCmp(Keypoint& f1, Keypoint& f2);
void Sift(const Mat &src, vector<Keypoint>& features, double sigma, int intervals=INTERVALS);

#endif