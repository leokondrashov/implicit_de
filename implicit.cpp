#include <iostream>
#include "vector.h"
#include "matrix.h"

#define tau 1e-3

double f1(const vector& u) {
	return 2 * u[0] * u[0] / (u[0] + 3) - 0.05 * u[0] - 5 * u[0] * u[1];
}

double f2(const vector& u) {
	return 0.0015 * u[0] * (1 - u[1] / 5) * (1 + (u[1] * u[1] / 0.0525 / 0.0525)) - 0.35 * u[1];
}

double J11(const vector& u) {
	return 2 * (u[0] * u[0] + 2 * u[0] * 3) / (u[0] + 3) / (u[0] + 3) - 0.05 - 5 * u[1];
}

double J12(const vector& u) {
	return -5 * u[0];
}

double J21(const vector& u) {
	return 0.0015 * (1 - u[1] / 5) * (1 + (u[1] * u[1] / 0.0525 / 0.0525));
}

double J22(const vector& u) {
	return -0.35 + 0.0015 * u[0] * (-1 / 5.0 + 2 * u[1] / 0.0525 / 0.0525 - 3 * u[1] * u[1] / 5 / 0.0525 / 0.0525);
}

double (*f[2])(const vector&) = {f1, f2};
double (*J[2][2])(const vector&) = {{J11, J12}, {J21, J22}};

vector RK1(const vector& u) { // explicit Euler method
	vector k(2);
	for (int i = 0; i < 3; i++)
		k[i] = f[i](u);
	return u + tau * k;
}

vector RK2(const vector& u) { // midpoint method 
	vector k[2] = {2, 2};

	for (int i = 0; i < 2; i++)
		k[0][i] = f[i](u);

	for (int i = 0; i < 2; i++)
		k[1][i] = f[i](u + tau / 2 * k[0]);

	return u + tau * k[1];
}

vector RK3(const vector& u) { // Heun's third order method
	vector k[3] = {2, 2, 2};
	matrix A(4); // Butcher tableau
	A[1][0] = 1 / 3.0; A[1][1] = 1 / 3.0;
	A[2][0] = 2 / 3.0; A[2][2] = 2 / 3.0;
	A[3][1] = 1 / 4.0; A[3][3] = 3 / 4.0; 

	for (int r = 0; r < 3; r++) {
		vector v = u;
		for (int i = 0 ; i < r; i++)
			v += tau * A[r][i + 1] * k[i];
		for (int i = 0; i < 2; i++)
			k[r][i] = f[i](v);
	}

	vector du(2);
	for (int i = 0; i < 3; i++)
		du += A[3][i + 1] * k[i];
	return u + tau * du;
}

vector RK4(const vector& u) { // Classic RK4
	vector k[4] = {2, 2, 2, 2};
	matrix A(5); // Butcher tableau
	A[1][0] = 1 / 2.0; A[1][1] = 1 / 2.0;
	A[2][0] = 1 / 2.0; A[2][2] = 1 / 2.0;
	A[3][0] = 1; A[3][3] = 1;
	A[4][1] = 1 / 6.0; A[4][2] = 1 / 3.0; A[4][3] = 1 / 3.0; A[4][4] = 1 / 6.0; 

	for (int r = 0; r < 4; r++) {
		vector v = u;
		for (int i = 0 ; i < r; i++)
			v += tau * A[r][i + 1] * k[i];
		for (int i = 0; i < 2; i++)
			k[r][i] = f[i](v);
	}

	vector du(2);
	for (int i = 0; i < 4; i++)
		du += A[4][i + 1] * k[i];
	return u + tau * du;
}

vector Rosenbrock(const vector& u) {
	matrix B(2);
	for (int i = 0; i < 2; i++) 
		for (int j = 0; j < 2; j++)
			B[i][j] = J[i][j](u);

	matrix E(2);
	E[0][0] = 1;
	E[1][1] = 1;

	vector F(2);
	for (int i = 0; i < 2; i++)
		F[i] = f[i](u);

	vector rhs(2);
	for (int i = 0; i < 2; i++)
		rhs[i] = f[i](u - 0.577 * tau * F);

	matrix lhs = E - 1.077 * tau * B + 0.372 * tau * tau * B * B;
	return u + tau * lhs.inverse() * rhs;
}

void dump(double t, const vector& u) {
	std::cout << t << ',' << u[0] << ',' << u[1] << std::endl;
}

int main() {
	double t = 0;
	vector u(2);
	u[0] = 0.8;
	u[1] = 0.05;
	//dump(t, u);

	for (double u00 = 2; u00 <= 4; u00 += 0.2) {
		for (double u10 = 0.15; u10 <= 0.25; u10 += 0.02) {
			t = 0;
			u[0] = u00;
			u[1] = u10;
			while (t <= 3) {
				u = Rosenbrock(u);
				t += tau;
				dump(t, u);
			}
		}
	}
/*
	while (t <= 50) {
		u = Rosenbrock(u);
		t += tau;
		dump(t, u);
	}
*/
}
