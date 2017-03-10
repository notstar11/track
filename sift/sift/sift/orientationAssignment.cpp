#include "orientationAssignment.h"

void CalculateScale(vector<Keypoint>& features, double sigma, int intervals)
{
	double intvl = 0;
	for (int i = 0; i < features.size(); i++)
	{
		intvl = features[i].interval + features[i].offset_interval;
		features[i].scale = sigma * pow(2.0, features[i].octave + intvl / intervals);
		features[i].octave_scale = sigma * pow(2.0, intvl / intervals);
	}

}

//对扩大的图像特征缩放
void HalfFeatures(vector<Keypoint>& features)
{
	for (int i = 0; i < features.size(); i++)
	{
		features[i].dx /= 2;
		features[i].dy /= 2;

		features[i].scale /= 2;
	}
}

bool CalcGradMagOri(const Mat& gauss, int x, int y, double& mag, double& ori)
{
	if (x > 0 && x < gauss.cols - 1 && y > 0 && y < gauss.rows - 1)
	{
		pixel_t *data = (pixel_t*)gauss.data;
		int step = gauss.step / sizeof(*data);

		double dx = *(data + step*y + (x + 1)) - (*(data + step*y + (x - 1)));
		double dy = *(data + step*(y + 1) + x) - (*(data + step*(y - 1) + x));

		mag = sqrt(dx*dx + dy*dy);

		//atan2返回[-Pi, -Pi]的弧度值
		ori = atan2(dy, dx);
		return true;
	}
	else
		return false;
}

double* CalculateOrientationHistogram(const Mat& gauss, int x, int y, int bins, int radius, double sigma)
{
	double *hist = new double[bins];

	for (int i = 0; i < bins; i++)
		*(hist + i) = 0.0;

	double mag, ori;

	double weight;

	int bin;
	const double PI2 = 2.0*CV_PI;

	double econs = -1.0 / (2.0*sigma*sigma);

	for (int i = -radius; i <= radius; i++)
	{
		for (int j = -radius; j <= radius; j++)
		{
			if (CalcGradMagOri(gauss, x + i, y + j, mag, ori))
			{
				weight = exp((i*i + j*j)*econs);

				//使用Pi-ori将ori转换到[0,2*PI]之间
				bin = cvRound(bins * (CV_PI - ori) / PI2);
				bin = bin < bins ? bin : 0;

				hist[bin] += mag * weight;
			}
		}
	}

	return hist;
}

//高斯平滑，模板为{0.25, 0.5, 0.25}
void GaussSmoothOriHist(double *hist, int n)
{
	double prev = hist[n - 1], temp, h0 = hist[0];


	for (int i = 0; i < n; i++)
	{
		temp = hist[i];
		hist[i] = 0.25 * prev + 0.5 * hist[i] +
			0.25 * (i + 1 >= n ? h0 : hist[i + 1]);
		prev = temp;
	}
}

//计算方向直方图中的主方向
double DominantDirection(double *hist, int n)
{
	double maxd = hist[0];
	for (int i = 1; i < n; i++)
	{
		if (hist[i] > maxd)
			maxd = hist[i];
	}
	return maxd;
}


void CopyKeypoint(const Keypoint& src, Keypoint& dst)
{
	dst.dx = src.dx;
	dst.dy = src.dy;

	dst.interval = src.interval;
	dst.octave = src.octave;
	dst.octave_scale = src.octave_scale;
	dst.offset_interval = src.offset_interval;

	dst.offset_x = src.offset_x;
	dst.offset_y = src.offset_y;

	dst.ori = src.ori;
	dst.scale = src.scale;
	dst.val = src.val;
	dst.x = src.x;
	dst.y = src.y;
}

//抛物插值


void CalcOriFeatures(const Keypoint& keypoint, vector<Keypoint>& features, const double *hist, int n, double mag_thr)
{
	double bin, PI2 = CV_PI * 2.0;
	int l, r;
	for (int i = 0; i < n; i++)
	{
		l = (i == 0) ? n - 1 : i - 1;
		r = (i + 1) % n;

		//hist[i]是极值
		if (hist[i] > hist[l] && hist[i] > hist[r] && hist[i] >= mag_thr)
		{
			bin = i + Parabola_Interpolate(hist[l], hist[i], hist[r]);
			bin = (bin < 0) ? (bin + n) : (bin >= n ? (bin - n) : bin);
			Keypoint new_key;
			CopyKeypoint(keypoint, new_key);
			new_key.ori = ((PI2 * bin) / n) - CV_PI;
			features.push_back(new_key);
		}
	}
}

//关键点方向分配
void OrientationAssignment(vector<Keypoint>& extrema, vector<Keypoint>& features, const vector<Mat>& gauss_pyr)
{
	int n = extrema.size();
	double *hist;

	for (int i = 0; i < n; i++)
	{

		hist = CalculateOrientationHistogram(gauss_pyr[extrema[i].octave*(INTERVALS + 3) + extrema[i].interval],
			extrema[i].x, extrema[i].y, ORI_HIST_BINS, cvRound(ORI_WINDOW_RADIUS*extrema[i].octave_scale),
			ORI_SIGMA_TIMES*extrema[i].octave_scale);

		for (int j = 0; j < ORI_SMOOTH_TIMES; j++)
			GaussSmoothOriHist(hist, ORI_HIST_BINS);
		double highest_peak = DominantDirection(hist, ORI_HIST_BINS);

		CalcOriFeatures(extrema[i], features, hist, ORI_HIST_BINS, highest_peak*ORI_PEAK_RATIO);

		delete[] hist;

	}
}