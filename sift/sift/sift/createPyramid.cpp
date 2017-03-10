#include "createPyramid.h"

//创建初始灰度图像
//初始图像先将原图像灰度化，再扩大一倍后，使用高斯模糊平滑
void CreateInitSmoothGray(const Mat &src, Mat &dst, double sigma)
{
	Mat gray, up;

	ConvertToGray(src, gray);
	//imshow("gray", gray);
	UpSample(gray, up);

	//-1层的sigma
	double  sigma_init = sqrt(sigma * sigma - (INIT_SIGMA * 2) * (INIT_SIGMA * 2));

	GaussianSmooth(up, dst, sigma_init);
}

//高斯金字塔
void GaussianPyramid(const Mat &src, vector<Mat>&gauss_pyr, int octaves, int intervals, double sigma)
{
	//
	double *sigmas = new double[intervals + 3];
	double k = pow(2.0, 1.0 / intervals);

	//cout <<"k=" <<k<<endl;
	sigmas[0] = sigma;
	/*
	for(int i = 1; i < intervals+3; i++)
	{
	sigmas[i] = k*sigmas[i-1];
	//cout << " "<<sigmas[i] ;
	}
	*/
	double sig_prev, sig_total;
	for (int i = 1; i < intervals + 3; i++)
	{
		sig_prev = pow(k, i - 1) * sigma;
		sig_total = sig_prev * k;
		sigmas[i] = sqrt(sig_total * sig_total - sig_prev * sig_prev);
	}

	for (int o = 0; o < octaves; o++)
	{
		//每组多三层
		for (int i = 0; i < intervals + 3; i++)
		{
			Mat mat;
			if (o == 0 && i == 0)
			{
				src.copyTo(mat);
			}
			else if (o != 0 && i == 0)
			{
				//前一组的倒数第二层
				DownSample(gauss_pyr[o*(intervals + 3) - 2], mat);
				//				DownSample(gauss_pyr[(o-1)*(intervals+3)+intervals], mat);
			}
			else
			{
				//每组中下一层由上一层高斯模糊得到
				GaussianSmooth(gauss_pyr[o * (intervals + 3) + i - 1], mat, sigmas[i]);
			}
			gauss_pyr.push_back(mat);
		}
	}

	delete[] sigmas;
}

//c = a - b
void Sub(const Mat& a, const Mat& b, Mat & c)
{
	if (a.rows != b.rows || a.cols != b.cols || a.type() != b.type())
		return;
	if (!c.empty())
		return;
	c.create(a.size(), a.type());

	pixel_t* ap = (pixel_t*)a.data;
	pixel_t* ap_end = (pixel_t*)a.dataend;
	pixel_t* bp = (pixel_t*)b.data;
	pixel_t* cp = (pixel_t*)c.data;
	int step = a.step / sizeof(pixel_t);

	while (ap != ap_end)
	{
		*cp++ = *ap++ - *bp++;
	}
}

//差分金字塔
void DogPyramid(const vector<Mat>& gauss_pyr, vector<Mat>& dog_pyr, int octaves, int intervals)
{
	for (int o = 0; o < octaves; o++)
	{
		for (int i = 1; i < intervals + 3; i++)
		{
			Mat mat;
			Sub(gauss_pyr[o*(intervals + 3) + i], gauss_pyr[o*(intervals + 3) + i - 1], mat);
			dog_pyr.push_back(mat);
		}
	}
}