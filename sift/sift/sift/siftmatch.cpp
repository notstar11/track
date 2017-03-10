#include "siftmatch.h"
void match(vector<Keypoint>&feature1,vector <Keypoint>&feature2,vector <pair<Keypoint,Keypoint>> &pair)
{
	int first = 0, second = 0;
	double firstMax =numeric_limits<double>::max() , secondMax = numeric_limits<double>::max();
	double value = 0.0;
	for (int i = 0; i < feature1.size(); i++){
		first = 0; second = 0;
		firstMax = numeric_limits<double>::max(); 
		secondMax = numeric_limits<double>::max();
		for (int j = 0; j < feature2.size(); j++){
			value = 0;
			for (int k = 0; k < FEATURE_ELEMENT_LENGTH; ++k)
			{
				value += (feature1[i].descriptor[k] - feature2[j].descriptor[k])*
					(feature1[i].descriptor[k] - feature2[j].descriptor[k]);
			}
			
			value = sqrt(value);
			if (value<=firstMax){
				firstMax = value;
				first = j;
			}
			else if (value <= secondMax){
				secondMax = value;
				second = j;
			}
			else{}
		}
		if (firstMax / secondMax < 0.48)
		{
			cout << (firstMax / secondMax) << endl;
			pair.push_back(make_pair(feature1[i], feature2[first]));
		}
			
	}
}
void showSift(const Mat &img1, const Mat &img2, Mat &dst, vector <pair<Keypoint, Keypoint>> &pair)
{
	img1.copyTo(dst(Rect(0, 0, img1.cols, img1.rows)));
	img2.copyTo(dst(Rect(img1.cols+200, 200, img2.cols, img2.rows)));
	for (int i = 0; i < pair.size(); ++i){
		Point dot1(pair[i].first.dx, pair[i].first.dy);
		Point dot2(pair[i].second.dx + img1.cols + 200, pair[i].second.dy+200);
		line(dst, dot1, dot2, Scalar(255, 0, 0));
	}
}