// [[Rcpp::depends(RcppEigen)]]
#include <iostream>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace std;
using namespace Eigen;

// [[Rcpp::export]]
double MV_Statistic(MatrixXd dist, int n, VectorXi group, int R)
{
	int *n_group = new int[R];
	for (int i = 0; i < R; i++)
	{
		n_group[i] = 0;
	}
	for (int i = 0; i < n; i++)
	{
		n_group[group(i) - 1]++;
	}
	double *F = new double[R + 1];
	double MV = 0;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			for (int r = 0; r < R + 1; r++)
			{
				F[r] = 0;
			}
			for (int k = 0; k < n; k++)
			{
				if (dist(i, k) <= dist(i, j))
				{
					F[0]++;
					F[group(k)]++;
				}
			}
			double tmp = 0;
			for (int r = 0; r < R; r++)
			{
				double ttmp = F[r + 1] / n_group[r] - F[0] / n;
				tmp += n_group[r] * ttmp * ttmp;
			}
			MV += tmp;
		}
	}
	MV = MV / n / n / n;
	return MV;
}