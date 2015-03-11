#define _CRT_SECURE_NO_WARNINGS

#include "calculation.h"
using namespace Garbage;

#include <ctime>
#include <iomanip>
#include <iostream>
using namespace std;

bool converged(matrix& x_prev, matrix& x_cur) {
	long double eps = 1e-3;
	long double norm = 0;
	for (int i = 0; i < x_prev.rows(); ++i) {
		norm += (x_prev[i][0] - x_cur[i][0])*(x_prev[i][0] - x_cur[i][0]);
	}
	return norm < eps;
}

int main() {
	freopen("output.out", "w", stdout);
	srand(time(nullptr));

	int n = 3;
	matrix b(n, 1);
	matrix x_prev(n, 1), x_cur(n, 1);
	matrix A(n, n);
	
	for (int i = 0; i < n; ++i) {
		b[i][0] = (i & 1) ? -1 : 1;
	}
	for (int i = 0; i < n; ++i) {
		x_prev[i][0] = 1;
	}
	for (int i = 0; i < n; ++i) {
		A[0][i] = 1;
	}
	for (int i = 1; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			A[i][j] = 1.0 / (i + j + 1);
		}
	}

	matrix LD(n, n), U(n, n);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j <= i; ++j) {
			LD[i][j] = A[i][j];
		}
	}
	for (int i = 0; i < n; ++i) {
		for (int j = i + 1; j < n; ++j) {
			U[i][j] = A[i][j];
		}
	}
	matrix inverse_LD = calculation::inverse_matrix::QR_method(LD);
	matrix T = inverse_LD * U;
	A = T * -1;

	int cnt = 0;
	int TL = 10000;

	while (cnt++ < TL) {
		/*for (int i = 0; i < n; ++i) {
			long double new_component = b[i][0];
			for (int j = 0; j < i; ++j) {
				new_component -= A[i][j] * x_cur[j][0];
			}
			for (int j = i + 1; j < n; ++j) {
				new_component -= A[i][j] * x_prev[j][0];
			}
			x_cur[i][0] = new_component / A[i][i];
		}*/
		x_cur = A * x_prev + inverse_LD * b;
		if (converged(x_prev, x_cur)) break;
		x_prev = x_cur;
		cout << x_cur << endl;
	}

	cout << fixed << setprecision(5) << x_cur << endl;
	return 0;
}