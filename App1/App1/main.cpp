#define _CRT_SECURE_NO_WARNINGS

#include "GaussSeidelMethod.h"

#include <algorithm>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <numeric>

using namespace std;
using namespace Garbage;

#define db(x) cout << #x << '=' << (x) << "\n"

void fill_matrix(matrix& A) {
	int n = A.rows();
	for (int i = 0; i < n; ++i) {
		A[0][i] = 1;
	}
	for (int i = 1; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			A[i][j] = long double(1.0) / long double(i + j + 1);
		}
	}
}

typedef long double ldouble;

int main() {
	freopen("output.out", "w", stdout);
	srand(time(nullptr));

	int n = 3;
	matrix A(n, n);

	fill_matrix(A);

	matrix b(n, 1);
	for (int i = 0; i < n; ++i) {
		b[i][0] = (i & 1) ? -1 : 1;
	}

	b = A.transpos() * b;
	A = A.transpos() * A;

	vector<ldouble> bb(n);
	for (int i = 0; i < n; ++i) {
		bb[i] = b[i][0];
	}

	int turns_limit = 10;
	int iter_cnt;
	int best_iter = 1e9;
	ldouble l = 0, r = 2;
	ldouble best_norm = 1e15, best_omega = -1;
	vector<ldouble> best_sol, x(n, 0), er, solution;
	for (int i = 1; i < turns_limit; ++i) {
		iter_cnt = 0;
		ldouble omega = l + (r - l) / turns_limit;
		db(omega);
		solution = GaussSeidelMethod::solve(A, x, bb, omega, iter_cnt);
		/*if (best_iter > iter_cnt) {
			best_iter = iter_cnt;
			best_omega = omega;
			best_sol = solution;
		}*/
		er = A * solution;
		ldouble norm = 0;
		for (int j = 0; j < n; ++j) {
			norm += (er[j] - bb[j]) * (er[j] - bb[j]);
		}
		if (best_norm > norm) {
			best_norm = norm;
			best_omega = omega;
			best_sol = solution;
		}
	}
	cout << iter_cnt << " " << best_omega << endl;
	for (int i = 0; i < n; ++i) {
		cout << best_sol[i] << " " << " ";
	}

	/*int iter_cnt = 0;
	ldouble omega = 1;
	vector<ldouble> x(n, 0);
	vector<ldouble> solution = GaussSeidelMethod::solve(A, x, bb, omega, iter_cnt);
	
	vector<ldouble> er = A * solution;
	cout << fixed << setprecision(10);
	for (int i = 0; i < n; ++i) {
		cout << er[i] << "\t" << bb[i] << endl;
	}*/

	return 0;
}