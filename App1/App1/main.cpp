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

int main() {
	freopen("output.out", "w", stdout);
	srand(time(nullptr));

	int n = 3;
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

	b = A.transpos() * b;
	A = A.transpos() * A;

	vector<long double> bb(n);
	for (int i = 0; i < n; ++i) {
		bb[i] = b[i][0];
	}

	int cnt = 10000;
	int anscnt = 10000;

	matrix best(n, 1);
	vector<long double> xi;
	
	long double minerr = 1e15;
	long double ansomega = -1;
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
			ansomega = omega;
			minerr = curerr;
			best = ans;
			anscnt = itercnt;
		}
	}

	cout << anscnt << endl;
	cout << fixed << setprecision(10) << ansomega << " " << minerr << endl;
	return 0;
}