#ifndef DISPLAYANDSAVE
#define DISPLAYANDSAVE
#include "basic.h"

void DrawSiftFeature(Mat& src, Keypoint& feat, CvScalar color);
void DrawSiftFeatures(Mat& src, vector<Keypoint>& features);
void DrawKeyPoints(Mat &src, vector<Keypoint>& keypoints);
const char* GetFileName(const char* dir, int i);
void cv64f_to_cv8U(const Mat& src, Mat& dst);
void writecv64f(const char* filename, const Mat& mat);
void write_pyr(const vector<Mat>& pyr, const char* dir);
void write_features(const vector<Keypoint>&features, const char*file);


#endif