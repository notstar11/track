#include "displayAndSave.h"

void DrawSiftFeature(Mat& src, Keypoint& feat, CvScalar color)
{
	int len, hlen, blen, start_x, start_y, end_x, end_y, h1_x, h1_y, h2_x, h2_y;
	double scl, ori;
	double scale = 5.0;
	double hscale = 0.75;
	CvPoint start, end, h1, h2;

	/* compute points for an arrow scaled and rotated by feat's scl and ori */
	start_x = cvRound(feat.dx);
	start_y = cvRound(feat.dy);
	scl = feat.scale;
	ori = feat.ori;
	len = cvRound(scl * scale);
	hlen = cvRound(scl * hscale);
	blen = len - hlen;
	end_x = cvRound(len *  cos(ori)) + start_x;
	end_y = cvRound(len * -sin(ori)) + start_y;
	h1_x = cvRound(blen *  cos(ori + CV_PI / 18.0)) + start_x;
	h1_y = cvRound(blen * -sin(ori + CV_PI / 18.0)) + start_y;
	h2_x = cvRound(blen *  cos(ori - CV_PI / 18.0)) + start_x;
	h2_y = cvRound(blen * -sin(ori - CV_PI / 18.0)) + start_y;
	start = cvPoint(start_x, start_y);
	end = cvPoint(end_x, end_y);
	h1 = cvPoint(h1_x, h1_y);
	h2 = cvPoint(h2_x, h2_y);

	line(src, start, end, color, 1, 8, 0);
	line(src, end, h1, color, 1, 8, 0);
	line(src, end, h2, color, 1, 8, 0);
}

void DrawSiftFeatures(Mat& src, vector<Keypoint>& features)
{
	CvScalar color = CV_RGB(0, 255, 0);
	for (int i = 0; i < features.size(); i++)
	{
		DrawSiftFeature(src, features[i], color);
	}
}

void DrawKeyPoints(Mat &src, vector<Keypoint>& keypoints)
{
	int j = 0;
	for (int i = 0; i < keypoints.size(); i++)
	{

		CvScalar color = { 255, 0, 0 };
		circle(src, Point(keypoints[i].dx, keypoints[i].dy), 2, color);
		j++;
	}
}

const char* GetFileName(const char* dir, int i)
{
	char *name = new char[50];
	sprintf(name, "%s\\%d\.jpg", dir, i);
	return name;
}

void cv64f_to_cv8U(const Mat& src, Mat& dst)
{
	double* data = (double*)src.data;
	int step = src.step / sizeof(*data);

	if (!dst.empty())
		return;
	dst.create(src.size(), CV_8U);

	uchar* dst_data = dst.data;

	for (int i = 0, m = 0; i < src.cols; i++, m++)
	{
		for (int j = 0, n = 0; j < src.rows; j++, n++)
		{
			*(dst_data + dst.step*j + i) = (uchar)(*(data + step*j + i) * 255);
		}
	}
}


//通过转换后保存的图像，会失真,和imshow显示出的图像相差很大
void writecv64f(const char* filename, const Mat& mat)
{
	Mat dst;
	cv64f_to_cv8U(mat, dst);
	imwrite(filename, dst);
}

void write_pyr(const vector<Mat>& pyr, const char* dir)
{
	for (int i = 0; i < pyr.size(); i++)
	{
		//		imshow(GetFileName("dog", i), dog_pyr[i]);
		//		imwrite(GetFileName("dogpyramid", i), dog_pyr[i]);

		writecv64f(GetFileName(dir, i), pyr[i]);
	}
}

void write_features(const vector<Keypoint>&features, const char*file)
{
	ofstream dout(file);
	dout << features.size() << endl;
	for (int i = 0; i < features.size(); i++)
	{
		dout << features[i].scale << " " << features[i].dx << " " << features[i].dy << endl;
		for (int j = 0; j < FEATURE_ELEMENT_LENGTH; j++)
		{
			if (j % 20 == 0)
				dout << endl;
			dout << features[i].descriptor[j] << " ";
		}
		dout << endl;
		dout << endl;
	}
}
