#ifndef _eigen_h_
#define _eigen_h_

#include "vector_matrix.h"

namespace nr
{
	extern "C" void zheev_(char *jobz, char *uplo, int *n, double *a, int *lda, double *w, double *work, int *lwork, double *rwork, int *info);
	
	inline void eig_sym(vec_doub &eigval, mat_complex &eigvec, const mat_complex &X)
	{
		char jobz = 'V';
		char uplo = 'U';
		
		if (X.nrows() != X.ncols()) throw("error: not a square matrix!");
		int n = X.nrows();
		
		double *a = new double[2*n*n];
		for (int j=0; j<n; j++)
		{
			for (int i=0; i<n; i++)			
			{
				a[2*(j*n+i)] = X(i,j).real();
				a[2*(j*n+i)+1] = X(i,j).imag();
			}
		}

		int lda = n;
		int lwork = 2*n-1;
		double *w = new double[n];
		double *work = new double[2*lwork];
		double *rwork = new double[3*n-2];
		int info = -1;

		zheev_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, rwork, &info);

		if (info == 0)
		{
			for (int i=0; i<n; i++)
			{
				eigval(i) = w[i];
				for (int j=0; j<n; j++)
				{
					eigvec(j,i).real(a[2*(i*n+j)]);
					eigvec(j,i).imag(a[2*(i*n+j)+1]);
				}
			}
		}
		else
		{
			throw("error: eigen value calculation failed!")
		}

		delete[] a, w, work, rwork;
		a = NULL;
		w = NULL;
		work = NULL;
		rwork = NULL; 
	}

}

#endif // _eigen_h_