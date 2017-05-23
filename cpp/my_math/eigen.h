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
		
		if (X.dim1() != X.dim2()) throw("error: not a square matrix!");
		int n = X.dim1();
		
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

	extern "C" void zheevx_(char *jobz, char *range, char *uplo, int *n, double *a, int *lda, double *vl, double *vu, int *il, int *iu, double *abstol, int *m, double *w, double *z, int *ldz, double *work, int *lwork, double *rwork, int *iwork, int *ifail, int *info);
	inline void eig_sym_selected(vec_doub &eigval, mat_complex &eigvec, const mat_complex &X, const int n_eigval)
	{
		char jobz = 'V';
		char range = 'I';
		char uplo = 'U';
		
		if (X.dim1() != X.dim2())
		{
			throw("error: not a square matrix!");
		}
		int n = X.dim1();
		
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
		double vl = 0;
		double vu = 1;
		int il = 1;
		int iu = n_eigval;
		double abstol = -1;
		int m = iu-il+1;
		double *w = new double[n];
		double *z = new double[2*n*m];
		int ldz = n;
		int lwork = 2*n+1; // this should be 2*n-1 but for some reason it does not work with that value and I have to use a larger value.
		double *work = new double[2*lwork];
		double *rwork = new double[7*n];
		int *iwork = new int[5*n];
		int *ifail = new int[n];
		int info = -1;

		zheevx_(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz, work, &lwork, rwork, iwork, ifail, &info);

		if (info == 0)
		{
			for (int i=0; i<m; i++)
			{
				eigval(i) = w[i];
				for (int j=0; j<n; j++)
				{
					eigvec(j,i).real(z[2*(i*n+j)]);
					eigvec(j,i).imag(z[2*(i*n+j)+1]);
				}
			}
		}
		else
		{
			throw("error: eigen value calculation failed!")
		}

		delete[] a, w, z, work, rwork, iwork, ifail;
		a = NULL;
		w = NULL;
		z = NULL;
		work = NULL;
		rwork = NULL;
		iwork = NULL;
		ifail = NULL;
	}

}

#endif // _eigen_h_