#define _CRT_SECURE_NO_WARNINGS

#include "calculation.h"
using namespace Garbage;

#define db(x) cout << #x << '=' << (x) << "\n"

#include <algorithm>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <numeric>
using namespace std;

class GaussSeidelMethod {

	typedef long double type;

public:

	static vector<type> solve(matrix& A, vector<type>& x, vector<type>& b, type omega, int& cntiter, int iteration_limit = 1000) {
		int n = x.size();

		type Cnorm;
		matrix DL(n, n), U(n, n);
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				if (j <= i) {
					DL[i][j] = A[i][j];
				}
				else {
					U[i][j] = A[i][j];
				}
			}
		}
		matrix C = calculation::inverse_matrix::LU_method(DL);
		C = C * U;
		Cnorm = norm(C);

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
			if (converged(A, Cnorm, x_prev, x_cur)) {
				db(it);
				break;
			}
			x_prev = x_cur;
			cntiter++;
		}
		return x_cur;
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
		type norm = 0;
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				norm += A[i][j] * A[i][j];
			}
		}
		return sqrt(norm);
	}

	static bool converged(matrix& A, vector<type>& x, vector<type>& b, type eps = 1e-3) {
		int n = A.rows();
		matrix xx(n, 1);
		for (int i = 0; i < n; ++i) {
			xx[i][0] = x[i];
		}
		matrix ans = A * xx;
		vector<type> err(n);
		for (int i = 0; i < n; ++i) {
			err[i] = ans[i][0] - b[i];
		}
		return norm(err) < eps;
	}

	static bool converged(matrix& A, type Cnorm, vector<type>& x_prev, vector<type>& x_cur, type eps = 1e-8) {
		int n = x_prev.size();
		vector<type> diff(n);
		for (int i = 0; i < n; ++i) {
			diff[i] = x_prev[i] - x_cur[i];
		}
		type errnorm = norm(diff);
		return abs((errnorm * Cnorm / (1 - Cnorm))) < eps;
	}
};

int main() {
	freopen("output.out", "w", stdout);
	srand(time(nullptr));

	int n = 100;
	matrix A(n, n);

	for (int i = 0; i < n; ++i) {
		A[0][i] = 1;
	}
	for (int i = 1; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			A[i][j] = long double(1.0) / long double(i + j + 1);
		}
	}

	matrix b(n, 1);
	vector<long double> x(n, 0);
	for (int i = 0; i < n; ++i) {
		b[i][0] = (i & 1) ? -1 : 1;
	}

	//matrix LUans = calculation::SLAE::LU_method(A, b);

	b = A.transpos() * b;
	A = A.transpos() * A;

	//cout << "A:" << endl << A << << endl;
	//cout << "b:" << endl << b << endl;

	vector<long double> bb(n);
	for (int i = 0; i < n; ++i) {
		bb[i] = b[i][0];
	}

	int anscnt;
	int cnt = 1000;
	vector<long double> xi;
	long double minerr = 1e15;
	matrix best(n, 1);
	long double l = 0, r = 2;
	for (int i = 1; i < cnt; ++i) {
		int itercnt = 0;
		long double omega = l + (r - l) * long double(i) / long double(cnt);
		xi = GaussSeidelMethod::solve(A, x, bb, omega, itercnt);
		matrix ans(n, 1);
		for (int j = 0; j < n; ++j) {
			ans[j][0] = xi[j];
		}
		matrix err = A * ans;
		long double curerr = 0;
		for (int j = 0; j < n; ++j) {
			curerr += (err[j][0] - bb[j]) * (err[j][0] - bb[j]);
		}
		if (minerr > curerr) {
			minerr = curerr;
			best = ans;
			anscnt = itercnt;
		}
	}

	cout << anscnt << endl;
	cout << fixed << setprecision(10) << minerr << endl;
	return 0;
}