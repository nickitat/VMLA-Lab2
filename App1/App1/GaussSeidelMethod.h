#ifndef __GAUSS_SEIDEL_METHOD__
#define __GAUSS_SEIDEL_METHOD__

#include "calculation.h"

#include <iostream>
#include <vector>

using namespace std;
using namespace Garbage;

#define db(x) cout << #x << '=' << (x) << "\n"

const int maxn = 100;
int C[maxn + 1][maxn + 1];

void C_calculate() {
	for (int n = 0; n <= maxn; ++n) {
		C[n][0] = C[n][n] = 1;
		for (int k = 1; k < n; ++k)
			C[n][k] = C[n - 1][k - 1] + C[n - 1][k];
	}
}

matrix get_inverse_matrix(int n) {
	C_calculate();
	matrix A(n, n);
	for (int i = 0; i < n; ++i) {
		A[i][0] = C[(n - 1) + i - 1][i - 1] * C[(n - 1)][i];
		if ((n - 1 - i) & 1) {
			A[i][0] *= -1;
		}
	}
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n - 1; ++j) {
			A[i][j + 1] = C[i + j][j] * C[i + j - 1][j - 1] * C[(n - 1) + i - 1][i + j] * C[(n - 1) + j][i + j] * i;
			if ((i - j) & 1) {
				A[i][j + 1] *= -1;
			}
		}
	}
	return A;
}

class GaussSeidelMethod {

	typedef long double type;

public:

	static vector<type> solve(matrix& A, vector<type>& x, vector<type>& b, type omega, int& cntiter, int iteration_limit = 1000000) {
		int n = x.size();

		matrix A_inv = get_inverse_matrix(n);
		type A_inv_norm = norm(A_inv);
		//db(A_inv_norm);

		vector<type> x_prev = x, x_cur(n);
		for (int it = 0; it < iteration_limit; ++it) {
			for (int i = 0; i < n; ++i) {
				type new_component = (1 - omega) * x_prev[i] + (omega)* b[i];
				for (int j = 0; j < i; ++j) {
					new_component -= A[i][j] * (omega)* x_cur[j];
				}
				for (int j = i + 1; j < n; ++j) {
					new_component -= A[i][j] * (omega)* x_prev[j];
				}
				x_cur[i] = new_component / A[i][i];
			}
			if (converged(A, A_inv_norm, x_cur, b)) {
				db(it);
				break;
			}
			x_prev = x_cur;
			cntiter++;
		}
		return vector<type>(x_cur);
	}

private:

	static type norm(vector<type>& v) {
		type norm = 0;
		for (int i = 0; i < v.size(); ++i) {
			norm += v[i] * v[i];
		}
		return sqrt(norm);
	}

	static type norm(matrix& A) {
		int n = A.rows();
		int m = A.columns();
		type norm = 0;
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) {
				norm += A[i][j] * A[i][j];
			}
		}
		return sqrt(norm);
	}

	static bool converged(matrix& A, type A_inv_norm, vector<type>& x_cur, vector<type>& b, type eps = 1e-8) {
		int n = x_cur.size();
		vector<type> r(n);
		for (int i = 0; i < n; ++i) {
			r[i] = -b[i];
			for (int j = 0; j < n; ++j) {
				r[i] += A[i][j] * x_cur[j];
			}
		}
		type nnorm = norm(r);
		return nnorm * A_inv_norm < eps;
	}
};

#endif/*__GAUSS_SEIDEL_METHOD__*/