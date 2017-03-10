//vs2010+opencv2.2
//zdd
//zddmail@gmail.com
//Just for fun


#include <windows.h>

#include <iostream>
using namespace std;

#include <cv.h>
#include <highgui.h>
#include <cxcore.h>
using namespace cv;

#include "sift.h"
#include "displayAndSave.h"
#include "siftmatch.h"


int main(int argc, char **argv)
{
	Mat src1 = imread("jobs_512.jpg");
	Mat src2 = imread("jobs.jpg");
	if(src1.empty()||src2.empty())
	{
		cout << "jobs_512.jpg or job.jpg open error! "<<endl;
		getchar();
		return -1;
	}

	if(src1.channels()!=3||src2.channels()!=3)
		return -2;

	vector<Keypoint> features1;
	vector<Keypoint>features2;

	Sift(src1, features1, 1.6);
	Sift(src2, features2, 1.6);
//	DrawKeyPoints(src, features);
	DrawSiftFeatures(src1, features1);
	DrawSiftFeatures(src2, features2);

	write_features(features1, "descriptor1.txt");
	write_features(features2, "descriptor1.txt");


	vector <pair<Keypoint, Keypoint>> pair;
	match(features1, features2, pair);
	Mat image(src1.rows+src2.rows+400,src1.cols+src2.cols+400,src1.type());
	showSift(src1, src2, image, pair);

	namedWindow("src1", WINDOW_AUTOSIZE);
	namedWindow("src2", WINDOW_AUTOSIZE);
	namedWindow("dst", WINDOW_AUTOSIZE);
	imshow("src1", src1);
	imshow("src2", src2);
	imshow("dst", image);
	waitKey(0);
	return 0;
}
