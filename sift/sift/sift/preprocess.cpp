#include "preprocess.h"
//转换为灰度图像
void ConvertToGray(const Mat& src, Mat& dst)
{
	Size size = src.size();
	if (dst.empty())
		dst.create(size, CV_64F);
	uchar* srcData = src.data;
	pixel_t* dstData = (pixel_t*)dst.data;
	int dstStep = dst.step / sizeof(dstData[0]);

	for (int j = 0; j < src.cols; j++)
	{
		for (int i = 0; i < src.rows; i++)
		{

			double b = *(srcData + src.step * i + src.channels() * j + 0) / 255.0;
			double g = *(srcData + src.step * i + src.channels() * j + 1) / 255.0;
			double r = *(srcData + src.step * i + src.channels() * j + 2) / 255.0;

			*((dstData + dstStep * i + dst.channels() * j)) = (b + g + r) / 3.0;
		}
	}
}

//隔点采样
void DownSample(const Mat& src, Mat& dst)
{
	if (src.channels() != 1)
		return;

	if (src.cols <= 1 || src.rows <= 1)
	{
		src.copyTo(dst);
		return;
	}

	dst.create((int)(src.rows / 2), (int)(src.cols / 2), src.type());
	//cout<<"-- "<<dst.rows<<" " <<dst.cols << " --"<<endl;

	pixel_t* srcData = (pixel_t*)src.data;
	pixel_t* dstData = (pixel_t*)dst.data;

	int srcStep = src.step / sizeof(srcData[0]);
	int dstStep = dst.step / sizeof(dstData[0]);

	int m = 0, n = 0;
	for (int j = 0; j < src.cols; j += 2, n++)
	{
		m = 0;
		for (int i = 0; i < src.rows; i += 2, m++)
		{
			pixel_t sample = *(srcData + srcStep * i + src.channels() * j);

			//防止当图像长宽不一致时，长宽为奇数时，m,n越界
			if (m < dst.rows && n < dst.cols)
			{
				*(dstData + dstStep * m + dst.channels() * n) = sample;

			}
		}
	}

}

//线性插值放大
void UpSample(const Mat &src, Mat &dst)
{
	if (src.channels() != 1)
		return;
	dst.create(src.rows * 2, src.cols * 2, src.type());

	pixel_t* srcData = (pixel_t*)src.data;
	pixel_t* dstData = (pixel_t*)dst.data;

	int srcStep = src.step / sizeof(srcData[0]);//每一维元素的个数
	int dstStep = dst.step / sizeof(dstData[0]);

	int m = 0, n = 0;
	for (int j = 0; j < src.cols - 1; j++, n += 2)
	{
		m = 0;
		for (int i = 0; i < src.rows - 1; i++, m += 2)
		{
			double sample = *(srcData + srcStep * i + src.channels() * j);
			*(dstData + dstStep * m + dst.channels() * n) = sample;

			double rs = *(srcData + srcStep * (i)+src.channels()*j) + (*(srcData + srcStep * (i + 1) + src.channels()*j));
			*(dstData + dstStep * (m + 1) + dst.channels() * n) = rs / 2;
			double cs = *(srcData + srcStep * i + src.channels()*(j)) + (*(srcData + srcStep * i + src.channels()*(j + 1)));
			*(dstData + dstStep * m + dst.channels() * (n + 1)) = cs / 2;

			double center = (*(srcData + srcStep * (i + 1) + src.channels() * j))
				+ (*(srcData + srcStep * i + src.channels() * j))
				+ (*(srcData + srcStep * (i + 1) + src.channels() * (j + 1)))
				+ (*(srcData + srcStep * i + src.channels() * (j + 1)));

			*(dstData + dstStep * (m + 1) + dst.channels() * (n + 1)) = center / 4;

		}

	}



	if (dst.rows < 3 || dst.cols < 3)
		return;

	//最后两行两列
	for (int k = dst.rows - 1; k >= 0; k--)
	{
		*(dstData + dstStep *(k)+dst.channels()*(dst.cols - 2)) = *(dstData + dstStep *(k)+dst.channels()*(dst.cols - 3));
		*(dstData + dstStep *(k)+dst.channels()*(dst.cols - 1)) = *(dstData + dstStep *(k)+dst.channels()*(dst.cols - 3));
	}
	for (int k = dst.cols - 1; k >= 0; k--)
	{
		*(dstData + dstStep *(dst.rows - 2) + dst.channels()*(k)) = *(dstData + dstStep *(dst.rows - 3) + dst.channels()*(k));
		*(dstData + dstStep *(dst.rows - 1) + dst.channels()*(k)) = *(dstData + dstStep *(dst.rows - 3) + dst.channels()*(k));
	}

}

//高斯平滑
//未使用sigma，边缘无处理
//未使用
void GaussianTemplateSmooth(const Mat &src, Mat &dst, double sigma)
{
	//高斯模板(7*7)，sigma = 0.84089642，归一化后得到
	static const double gaussianTemplate[7][7] =
	{
		{ 0.00000067, 0.00002292, 0.00019117, 0.00038771, 0.00019117, 0.00002292, 0.00000067 },
		{ 0.00002292, 0.00078633, 0.00655965, 0.01330373, 0.00655965, 0.00078633, 0.00002292 },
		{ 0.00019117, 0.00655965, 0.05472157, 0.11098164, 0.05472157, 0.00655965, 0.00019117 },
		{ 0.00038771, 0.01330373, 0.11098164, 0.22508352, 0.11098164, 0.01330373, 0.00038771 },
		{ 0.00019117, 0.00655965, 0.05472157, 0.11098164, 0.05472157, 0.00655965, 0.00019117 },
		{ 0.00002292, 0.00078633, 0.00655965, 0.01330373, 0.00655965, 0.00078633, 0.00002292 },
		{ 0.00000067, 0.00002292, 0.00019117, 0.00038771, 0.00019117, 0.00002292, 0.00000067 }
	};

	dst.create(src.size(), src.type());
	uchar* srcData = src.data;
	uchar* dstData = dst.data;

	for (int j = 0; j < src.cols - 7; j++)
	{
		for (int i = 0; i < src.rows - 7; i++)
		{
			double acc = 0;
			double accb = 0, accg = 0, accr = 0;
			for (int m = 0; m < 7; m++)
			{
				for (int n = 0; n < 7; n++)
				{
					if (src.channels() == 1)
						acc += *(srcData + src.step * (i + n) + src.channels() * (j + m)) * gaussianTemplate[m][n];
					else
					{
						accb += *(srcData + src.step * (i + n) + src.channels() * (j + m) + 0) * gaussianTemplate[m][n];
						accg += *(srcData + src.step * (i + n) + src.channels() * (j + m) + 1) * gaussianTemplate[m][n];
						accr += *(srcData + src.step * (i + n) + src.channels() * (j + m) + 2) * gaussianTemplate[m][n];
					}
				}
			}
			if (src.channels() == 1)
				*(dstData + dst.step * (i + 3) + dst.channels() * (j + 3)) = (int)acc;
			else
			{
				*(dstData + dst.step * (i + 3) + dst.channels() * (j + 3) + 0) = (int)accb;
				*(dstData + dst.step * (i + 3) + dst.channels() * (j + 3) + 1) = (int)accg;
				*(dstData + dst.step * (i + 3) + dst.channels() * (j + 3) + 2) = (int)accr;
			}
		}
	}

}

//未使用
void GaussianSmooth2D(const Mat &src, Mat &dst, double sigma)
{
	if (src.channels() != 1)
		return;

	//确保sigma为正数 
	sigma = sigma > 0 ? sigma : 0;
	//高斯核矩阵的大小为(6*sigma+1)*(6*sigma+1)
	//ksize为奇数
	int ksize = cvRound(sigma * 3) * 2 + 1;

	if (ksize == 1)
	{
		src.copyTo(dst);
		return;
	}

	dst.create(src.size(), src.type());

	//计算高斯核矩阵
	double *kernel = new double[ksize*ksize];

	double scale = -0.5 / (sigma*sigma);
	const double PI = 3.141592653;
	double cons = -scale / PI;

	double sum = 0;

	for (int i = 0; i < ksize; i++)
	{
		for (int j = 0; j < ksize; j++)
		{
			int x = i - (ksize - 1) / 2;
			int y = j - (ksize - 1) / 2;
			kernel[i*ksize + j] = cons * exp(scale * (x*x + y*y));
			sum += kernel[i*ksize + j];
		}
	}
	//归一化
	for (int i = ksize*ksize - 1; i >= 0; i--)
	{
		*(kernel + i) /= sum;
	}
	uchar* srcData = src.data;
	uchar* dstData = dst.data;

	//图像卷积运算
	for (int j = 0; j < src.cols - ksize; j++)
	{
		for (int i = 0; i < src.rows - ksize; i++)
		{
			double acc = 0;

			for (int m = 0; m < ksize; m++)
			{
				for (int n = 0; n < ksize; n++)
				{
					acc += *(srcData + src.step * (i + n) + src.channels() * (j + m)) * kernel[m*ksize + n];
				}
			}

			/*
			for(int l = 0; l < ksize * ksize; l++)
			acc +=  *(srcData + src.step * (i+(int)l/ksize) + src.channels() * (j+(int)l%ksize)) * kernel[l];
			*/
			*(dstData + dst.step * (i + (ksize - 1) / 2) + (j + (ksize - 1) / 2)) = (int)acc;
		}
	}
	delete[]kernel;
}

void GaussianSmooth(const Mat &src, Mat &dst, double sigma)
{
	GaussianBlur(src, dst, Size(0, 0), sigma);
}