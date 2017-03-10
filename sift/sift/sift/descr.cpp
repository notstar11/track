#include "descr.h"
void InterpHistEntry(double ***hist, double xbin, double ybin, double obin, double mag, int bins, int d)
{
	double d_r, d_c, d_o, v_r, v_c, v_o;
	double** row, *h;
	int r0, c0, o0, rb, cb, ob, r, c, o;

	r0 = cvFloor(ybin);
	c0 = cvFloor(xbin);
	o0 = cvFloor(obin);
	d_r = ybin - r0;
	d_c = xbin - c0;
	d_o = obin - o0;

	/*
	做插值：
	xbin,ybin,obin:种子点所在子窗口的位置和方向
	所有种子点都将落在4*4的窗口中
	r0,c0取不大于xbin，ybin的正整数
	r0,c0只能取到0,1,2
	xbin,ybin在(-1, 2)

	r0取不大于xbin的正整数时。
	r0+0 <= xbin <= r0+1
	mag在区间[r0,r1]上做插值

	obin同理
	*/

	for (r = 0; r <= 1; r++)
	{
		rb = r0 + r;
		if (rb >= 0 && rb < d)
		{
			v_r = mag * ((r == 0) ? 1.0 - d_r : d_r);
			row = hist[rb];
			for (c = 0; c <= 1; c++)
			{
				cb = c0 + c;
				if (cb >= 0 && cb < d)
				{
					v_c = v_r * ((c == 0) ? 1.0 - d_c : d_c);
					h = row[cb];
					for (o = 0; o <= 1; o++)
					{
						ob = (o0 + o) % bins;
						v_o = v_c * ((o == 0) ? 1.0 - d_o : d_o);
						h[ob] += v_o;
					}
				}
			}
		}
	}


}

double*** CalculateDescrHist(const Mat& gauss, int x, int y, double octave_scale, double ori, int bins, int width)
{
	double ***hist = new double**[width];

	//申请空间并初始化
	for (int i = 0; i < width; i++)
	{
		hist[i] = new double*[width];
		for (int j = 0; j < width; j++)
		{
			hist[i][j] = new double[bins];
		}
	}

	for (int r = 0; r < width; r++)
	for (int c = 0; c < width; c++)
	for (int o = 0; o < bins; o++)
		hist[r][c][o] = 0.0;


	double cos_ori = cos(ori);
	double sin_ori = sin(ori);

	//6.1高斯权值，sigma等于描述字窗口宽度的一半
	double sigma = 0.5 * width;
	double conste = -1.0 / (2 * sigma*sigma);

	double PI2 = CV_PI * 2;

	double sub_hist_width = DESCR_SCALE_ADJUST * octave_scale;

	//领域半径
	int radius = (sub_hist_width*sqrt(2.0)*(width + 1)) / 2.0 + 0.5; //+0.5取四舍五入

	double grad_ori, grad_mag;
	for (int i = -radius; i <= radius; i++)
	{
		for (int j = -radius; j <= radius; j++)
		{
			double rot_x = (cos_ori * j - sin_ori * i) / sub_hist_width;
			double rot_y = (sin_ori * j + cos_ori * i) / sub_hist_width;


			//xbin,ybin为落在4*4窗口中的下标值
			double xbin = rot_x + width / 2 - 0.5;
			double ybin = rot_y + width / 2 - 0.5;

			//
			if (xbin > -1.0 && xbin < width && ybin > -1.0 && ybin < width)
			{
				if (CalcGradMagOri(gauss, x + j, y + i, grad_mag, grad_ori))
				{
					grad_ori = (CV_PI - grad_ori) - ori;
					while (grad_ori < 0.0)
						grad_ori += PI2;
					while (grad_ori >= PI2)
						grad_ori -= PI2;

					double obin = grad_ori * (bins / PI2);

					double weight = exp(conste*(rot_x*rot_x + rot_y * rot_y));

					InterpHistEntry(hist, xbin, ybin, obin, grad_mag*weight, bins, width);

				}
			}
		}
	}

	return hist;
}

void NormalizeDescr(Keypoint& feat)
{
	double cur, len_inv, len_sq = 0.0;
	int i, d = feat.descr_length;

	for (i = 0; i < d; i++)
	{
		cur = feat.descriptor[i];
		len_sq += cur*cur;
	}
	len_inv = 1.0 / sqrt(len_sq);
	for (i = 0; i < d; i++)
		feat.descriptor[i] *= len_inv;
}

void HistToDescriptor(double ***hist, int width, int bins, Keypoint& feature)
{
	int int_val, i, r, c, o, k = 0;

	for (r = 0; r < width; r++)
	for (c = 0; c < width; c++)
	for (o = 0; o < bins; o++)
	{
		feature.descriptor[k++] = hist[r][c][o];
	}

	feature.descr_length = k;
	NormalizeDescr(feature);
	for (i = 0; i < k; i++)
	if (feature.descriptor[i] > DESCR_MAG_THR)
		feature.descriptor[i] = DESCR_MAG_THR;
	NormalizeDescr(feature);

	/* convert floating-point descriptor to integer valued descriptor */
	for (i = 0; i < k; i++)
	{
		int_val = INT_DESCR_FCTR * feature.descriptor[i];
		feature.descriptor[i] = min(255, int_val);
	}
}

//计算描述符
void DescriptorRepresentation(vector<Keypoint>& features, const vector<Mat>& gauss_pyr, int bins, int width)
{
	double ***hist;

	for (int i = 0; i < features.size(); i++)
	{
		hist = CalculateDescrHist(gauss_pyr[features[i].octave*(INTERVALS + 3) + features[i].interval],
			features[i].x, features[i].y, features[i].octave_scale, features[i].ori, bins, width);


		HistToDescriptor(hist, width, bins, features[i]);

		for (int j = 0; j < width; j++)
		{

			for (int k = 0; k < width; k++)
			{
				delete[] hist[j][k];
			}
			delete[] hist[j];
		}
		delete[] hist;
	}
}
