#include "sift.h"

bool FeatureCmp(Keypoint& f1, Keypoint& f2)
{
	return f1.scale < f2.scale;
}

//sift 算法
void Sift(const Mat &src, vector<Keypoint>& features, double sigma, int intervals)
{
	Mat init_gray;
	CreateInitSmoothGray(src, init_gray, sigma);
	//计算组数
	//int octaves = log((double)min(init_gray.rows, init_gray.cols))/log(2.0) - 3;
	int octaves = log((double)min(init_gray.rows, init_gray.cols)) / log(2.0) - 2;
	cout << "rows=" << init_gray.rows << " cols=" << init_gray.cols << " octaves=" << octaves << endl;


	cout << "building gaussian pyramid ..." << endl;
	vector<Mat> gauss_pyr;
	GaussianPyramid(init_gray, gauss_pyr, octaves, intervals, sigma);

	//	write_pyr(gauss_pyr, "gausspyramid");

	cout << "building difference of gaussian pyramid..." << endl;
	vector<Mat> dog_pyr;
	DogPyramid(gauss_pyr, dog_pyr, octaves, intervals);

	//	write_pyr(dog_pyr, "dogpyramid");

	cout << "deatecting local extrema..." << endl;

	vector<Keypoint> extrema;
	DetectionLocalExtrema(dog_pyr, extrema, octaves, intervals);

	cout << "--keypoints cout: " << extrema.size() << " --" << endl;

	cout << "extrema detection finished." << endl << "--please look dir gausspyramid, dogpyramid and extrema.txt.--" << endl;
	//计算尺度
	CalculateScale(extrema, sigma, intervals);

	HalfFeatures(extrema);

	cout << "orientation assignment..." << endl;
	OrientationAssignment(extrema, features, gauss_pyr);
	cout << "--features count: " << features.size() << " --" << endl;

	DescriptorRepresentation(features, gauss_pyr, DESCR_HIST_BINS, DESCR_WINDOW_WIDTH);

	sort(features.begin(), features.end(), FeatureCmp);

	cout << "finished." << endl;
}
