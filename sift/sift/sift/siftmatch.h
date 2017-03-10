#ifndef SIFTMATCH
#define SIFTMATCH
#include "basic.h"

void match(vector<Keypoint>&feature1, vector <Keypoint>&feature2, vector <pair<Keypoint, Keypoint>> &pair);
void showSift(const Mat &img1, const Mat &img2, Mat &dst, vector <pair<Keypoint, Keypoint>> &pair);


#endif