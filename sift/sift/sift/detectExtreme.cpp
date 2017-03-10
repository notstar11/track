#include "detectExtreme.h"
//
bool isExtremum(int x, int y, const vector<Mat>& dog_pyr, int index)
{
	pixel_t * data = (pixel_t *)dog_pyr[index].data;
	int step = dog_pyr[index].step / sizeof(data[0]);
	pixel_t val = *(data + y*step + x);

	if (val > 0)
	{
		for (int i = -1; i <= 1; i++)
		{
			int stp = dog_pyr[index + i].step / sizeof(pixel_t);
			for (int j = -1; j <= 1; j++)
			{
				for (int k = -1; k <= 1; k++)
				{
					//检查最大极值
					if (val < *((pixel_t*)dog_pyr[index + i].data + stp*(y + j) + (x + k)))
					{
						return false;
					}
				}
			}
		}
	}
	else
	{
		for (int i = -1; i <= 1; i++)
		{
			int stp = dog_pyr[index + i].step / sizeof(pixel_t);
			for (int j = -1; j <= 1; j++)
			{
				for (int k = -1; k <= 1; k++)
				{
					//检查最小极值
					if (val > *((pixel_t*)dog_pyr[index + i].data + stp*(y + j) + (x + k)))
					{
						return false;
					}
				}
			}
		}
	}

	return true;
}

//4.1 eliminating edge responses
//hessian矩阵，排除边缘点 
bool passEdgeResponse(int x, int y, const vector<Mat>& dog_pyr, int index, double r)
{
	pixel_t *data = (pixel_t *)dog_pyr[index].data;
	int step = dog_pyr[index].step / sizeof(data[0]);
	pixel_t val = *(data + y*step + x);

	double Dxx, Dyy, Dxy;
	double Tr_h, Det_h;

	//hessian矩阵
	//	   _ 	    _
	//    | Dxx  Dxy |
	// H =|			 |
	//	  |_Dxy  Dyy_|	
	//	  
	Dxx = DAt(x + 1, y) + DAt(x - 1, y) - 2 * val;
	Dyy = DAt(x, y + 1) + DAt(x, y - 1) - 2 * val;
	Dxy = (DAt(x + 1, y + 1) + DAt(x - 1, y - 1) - DAt(x - 1, y + 1) - DAt(x + 1, y - 1)) / 4.0;

	Tr_h = Dxx + Dyy;
	Det_h = Dxx * Dyy - Dxy * Dxy;

	if (Det_h <= 0)
		return false;

	if (Tr_h * Tr_h / Det_h < (r + 1) * (r + 1) / r)
		return true;

	return false;
}


//#define At(index, x, y) (*((pixel_t*)dog_pyr[(index)].data+(y)*((int)(dog_pyr[(index)].step/sizeof((pixel_t*)dog_pyr[index].data[0])))+(x)))
double PyrAt(const vector<Mat>& pyr, int index, int x, int y)
{
	pixel_t *data = (pixel_t*)pyr[index].data;
	int step = pyr[index].step / sizeof(data[0]);
	pixel_t val = *(data + y*step + x);

	return val;
}


//3维D(x)一阶偏导,dx列向量
void DerivativeOf3D(int x, int y, const vector<Mat>& dog_pyr, int index, double *dx)
{
	double Dx = (At(index, x + 1, y) - At(index, x - 1, y)) / 2.0;
	double Dy = (At(index, x, y + 1) - At(index, x, y - 1)) / 2.0;
	double Ds = (At(index + 1, x, y) - At(index - 1, x, y)) / 2.0;

	dx[0] = Dx;
	dx[1] = Dy;
	dx[2] = Ds;
}

//3维D(x)二阶偏导，即Hessian矩阵
void Hessian3D(int x, int y, const vector<Mat>& dog_pyr, int index, double *H)
{
	double val, Dxx, Dyy, Dss, Dxy, Dxs, Dys;

	val = At(index, x, y);

	Dxx = At(index, x + 1, y) + At(index, x - 1, y) - 2 * val;
	Dyy = At(index, x, y + 1) + At(index, x, y - 1) - 2 * val;
	Dss = At(index + 1, x, y) + At(index - 1, x, y) - 2 * val;

	Dxy = (At(index, x + 1, y + 1) + At(index, x - 1, y - 1)
		- At(index, x + 1, y - 1) - At(index, x - 1, y + 1)) / 4.0;

	Dxs = (At(index + 1, x + 1, y) + At(index - 1, x - 1, y)
		- At(index - 1, x + 1, y) - At(index + 1, x - 1, y)) / 4.0;

	Dys = (At(index + 1, x, y + 1) + At(index - 1, x, y - 1)
		- At(index + 1, x, y - 1) - At(index - 1, x, y + 1)) / 4.0;

	Hat(0, 0) = Dxx;
	Hat(1, 1) = Dyy;
	Hat(2, 2) = Dss;

	Hat(1, 0) = Hat(0, 1) = Dxy;
	Hat(2, 0) = Hat(0, 2) = Dxs;
	Hat(2, 1) = Hat(1, 2) = Dys;
}

//3*3阶矩阵求逆
bool Inverse3D(const double *H, double *H_inve)
{
	//A=|H|
	//		 / A00 A01 A02 \				   
	//若H =  | A10 A11 A12 |   
	//		 \ A20 A21 A22 /	
	//则 行列式|H|=A00*A11*A22+A01*A12*A20+A02*A10*A21
	//	    -A00*A12*A21-A01*A10*A22-A02*A11*A20
	//

	double A = Hat(0, 0)*Hat(1, 1)*Hat(2, 2)
		+ Hat(0, 1)*Hat(1, 2)*Hat(2, 0)
		+ Hat(0, 2)*Hat(1, 0)*Hat(2, 1)
		- Hat(0, 0)*Hat(1, 2)*Hat(2, 1)
		- Hat(0, 1)*Hat(1, 0)*Hat(2, 2)
		- Hat(0, 2)*Hat(1, 1)*Hat(2, 0);
	//cout<<A<<endl;
	//没有逆矩阵
	if (fabs(A) < 1e-10)
		return false;



	//三阶逆矩阵运算公式：
	//		 / a b c \				    / ei-hf -(bi-ch) bf-ce\
		//若A =  | d e f |   则A(-1) =1/|H|*| fg-id -(cg-ia) cd-af |
	//		 \ g h i /				    \ dh-ge -(ah-gb) ae-bd/



	HIat(0, 0) = Hat(1, 1) * Hat(2, 2) - Hat(2, 1)*Hat(1, 2);
	HIat(0, 1) = -(Hat(0, 1) * Hat(2, 2) - Hat(2, 1) * Hat(0, 2));
	HIat(0, 2) = Hat(0, 1) * Hat(1, 2) - Hat(0, 2)*Hat(1, 1);

	HIat(1, 0) = Hat(1, 2) * Hat(2, 0) - Hat(2, 2)*Hat(1, 0);
	HIat(1, 1) = -(Hat(0, 2) * Hat(2, 0) - Hat(0, 0) * Hat(2, 2));
	HIat(1, 2) = Hat(0, 2) * Hat(1, 0) - Hat(0, 0)*Hat(1, 2);

	HIat(2, 0) = Hat(1, 0) * Hat(2, 1) - Hat(1, 1)*Hat(2, 0);
	HIat(2, 1) = -(Hat(0, 0) * Hat(2, 1) - Hat(0, 1) * Hat(2, 0));
	HIat(2, 2) = Hat(0, 0) * Hat(1, 1) - Hat(0, 1)*Hat(1, 0);

	for (int i = 0; i < 9; i++)
	{
		*(H_inve + i) /= A;
	}
	return true;
}

//计算x^
void GetOffsetX(int x, int y, const vector<Mat>& dog_pyr, int index, double *offset_x)
{
	//x^ = -H^(-1) * dx; dx = (Dx, Dy, Ds)^T
	double H[9], H_inve[9] = { 0 };
	Hessian3D(x, y, dog_pyr, index, H);
	Inverse3D(H, H_inve);
	double dx[3];
	DerivativeOf3D(x, y, dog_pyr, index, dx);

	for (int i = 0; i < 3; i++)
	{
		offset_x[i] = 0.0;
		for (int j = 0; j < 3; j++)
		{
			offset_x[i] += H_inve[i * 3 + j] * dx[j];
		}
		offset_x[i] = -offset_x[i];
	}
}

//计算|D(x^)|
double GetFabsDx(int x, int y, const vector<Mat>& dog_pyr, int index, const double* offset_x)
{
	//|D(x^)|=D + 0.5 * dx * offset_x; dx=(Dx, Dy, Ds)^T
	double dx[3];
	DerivativeOf3D(x, y, dog_pyr, index, dx);

	double term = 0.0;
	for (int i = 0; i < 3; i++)
		term += dx[i] * offset_x[i];

	pixel_t *data = (pixel_t *)dog_pyr[index].data;
	int step = dog_pyr[index].step / sizeof(data[0]);
	pixel_t val = *(data + y*step + x);

	return fabs(val + 0.5 * term);
}

//修正极值点，删除不稳定点
// |D(x)| < 0.03 Lowe 2004
Keypoint* InterploationExtremum(int x, int y, const vector<Mat>& dog_pyr, int index, int octave, int interval, double dxthreshold)
{
	//计算x=(x,y,sigma)^T
	//x^ = -H^(-1) * dx; dx = (Dx, Dy, Ds)^T
	double offset_x[3] = { 0 };
	//
	const Mat &mat = dog_pyr[index];
	int idx = index;
	int intvl = interval;
	int i = 0;
	while (i < MAX_INTERPOLATION_STEPS)
	{
		GetOffsetX(x, y, dog_pyr, idx, offset_x);
		//4. Accurate keypoint localization.  Lowe
		//
		//如果offset_x 的任一维度大于0.5，it means that the extremum lies closer to a different sample point.
		if (fabs(offset_x[0]) < 0.5 && fabs(offset_x[1]) < 0.5 && fabs(offset_x[2]) < 0.5)
			break;

		//用周围的点代替
		//
		x += cvRound(offset_x[0]);
		y += cvRound(offset_x[1]);
		interval += cvRound(offset_x[2]);

		idx = index - intvl + interval;
		//		idx = octave*(INTERVALS+2)+interval;

		if (interval < 1 || interval > INTERVALS ||
			x >= mat.cols - 1 || x < 2 ||
			y >= mat.rows - 1 || y < 2)  //此处保证检测边时 x+1,y+1和x-1, y-1有效
		{
			return NULL;
		}

		i++;
	}

	//窜改失败
	if (i >= MAX_INTERPOLATION_STEPS)
		return NULL;

	//rejecting unstable extrema
	//|D(x^)| < 0.03取经验值
	if (GetFabsDx(x, y, dog_pyr, idx, offset_x) < dxthreshold / INTERVALS)
	{
		return NULL;
	}

	Keypoint *keypoint = new Keypoint;


	keypoint->x = x;
	keypoint->y = y;

	keypoint->offset_x = offset_x[0];
	keypoint->offset_y = offset_x[1];

	keypoint->interval = interval;
	keypoint->offset_interval = offset_x[2];

	keypoint->octave = octave;

	keypoint->dx = (x + offset_x[0])*pow(2.0, octave);
	keypoint->dy = (y + offset_x[1])*pow(2.0, octave);

	return keypoint;
}

//检测当地极值点
void DetectionLocalExtrema(const vector<Mat>& dog_pyr, vector<Keypoint>& extrema, int octaves, int intervals)
{
	long int dd = 0, cc1 = 0, cc2 = 0, cc3 = 0, cc0 = 0, cc00 = 0;

	double thresh = 0.5 * DXTHRESHOLD / intervals;
	for (int o = 0; o < octaves; o++)
	{
		//第一层和最后一层极值忽略
		for (int i = 1; i < (intervals + 2) - 1; i++)
		{
			int index = o*(intervals + 2) + i;
			pixel_t *data = (pixel_t *)dog_pyr[index].data;
			int step = dog_pyr[index].step / sizeof(data[0]);


			for (int y = 1; y < dog_pyr[index].rows - 2; y++)
			{
				for (int x = 1; x < dog_pyr[index].cols - 2; x++)
				{
					cc00++;
					//
					pixel_t val = *(data + y*step + x);
					if (fabs(val) > thresh) //排除阈值过小的点
					{
						cc0++;
						if (isExtremum(x, y, dog_pyr, index))
						{
							cc1++;
							Keypoint *extrmum = InterploationExtremum(x, y, dog_pyr, index, o, i);

							if (extrmum)
							{
								cc2++;

								if (passEdgeResponse(extrmum->x, extrmum->y, dog_pyr, index))
								{
									extrmum->val = *(data + extrmum->y*step + extrmum->x);
									cc3++;
									extrema.push_back(*extrmum);
								}

								delete extrmum;
							}
						}
					}
				}
			}

		}
	}
	cout << "-- " << "cc00: " << cc00 << ", cc0: " << cc0 << ", cc1: " << cc1 << ", cc2: " << cc2 << ", cc3: " << cc3 << " " << thresh << " --" << endl;
}
