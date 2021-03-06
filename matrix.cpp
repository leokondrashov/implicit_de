#include <iostream>
#include <cmath>

#include "matrix.h"

matrix::matrix(int n) : n(n) {
	lines.resize(n, {n});
}

matrix::matrix(const matrix& oth) {
	n = oth.n;
	lines = oth.lines;
}

matrix::matrix(matrix&& oth) {
	n = oth.n;
	oth.n = 0;
	lines = std::move(oth.lines);
}

matrix& matrix::operator=(const matrix& oth) {
	n = oth.n;
	lines = oth.lines;

	return *this;
}

matrix& matrix::operator=(matrix&& oth) {
	n = oth.n;
	oth.n = 0;
	lines = std::move(oth.lines);

	return *this;
}

matrix::~matrix() {
}

matrix::matrix_line::matrix_line(int n) : n(n){
	data.resize(n);
}

matrix::matrix_line::matrix_line(const matrix_line& ml) {
	n = ml.n;
	data = ml.data;
}

matrix::matrix_line::matrix_line(matrix_line&& ml) {
	n = ml.n;
	ml.n = 0;
	data = std::move(ml.data);
}

matrix::matrix_line& matrix::matrix_line::operator=(const matrix_line& ml) {
	n = ml.n;
	data = ml.data;

	return *this;
}

matrix::matrix_line& matrix::matrix_line::operator=(matrix_line&& ml) {
	n = ml.n;
	ml.n = 0;
	data = std::move(ml.data);

	return *this;
}

matrix::matrix_line::~matrix_line() {
}

const matrix::matrix_line& matrix::operator[](int i) const {
	return lines[i];
}

const double& matrix::matrix_line::operator[](int i) const {
	return data[i];
}

matrix::matrix_line& matrix::operator[](int i) {
	return lines[i];
}

double& matrix::matrix_line::operator[](int i) {
	return data[i];
}

void matrix::dump() {
	std::cout << "matrix [" << this << "]:" << std::endl;
	for (int i = 0; i < n; i++) {
		std::cout << "[ ";
		for (int j = 0; j < n; j++) {
			std::cout << (*this)[i][j] << " ";
		}
		std::cout << "]" << std::endl;
	}
}

matrix::matrix_line& matrix::matrix_line::operator+=(const matrix_line& ml) {
	for (int i = 0; i < n; i++) {
		data[i] += ml.data[i];
	}

	return *this;
}

matrix::matrix_line& matrix::matrix_line::operator*=(double d) {
	for (int i = 0; i < n; i++) {
		data[i] *= d;
	}

	return *this;
}

matrix::matrix_line& matrix::matrix_line::operator/=(double d) {
	for (int i = 0; i < n; i++) {
		data[i] /= d;
	}

	return *this;
}

matrix::matrix_line operator+(const matrix::matrix_line& left, const matrix::matrix_line& right) {
	matrix::matrix_line res{left};
	return res += right;
}

matrix::matrix_line operator*(const matrix::matrix_line& left, double right) {
	matrix::matrix_line res{left};
	return res *= right;
}

matrix::matrix_line operator*(double left, const matrix::matrix_line& right) {
	matrix::matrix_line res{right};
	return res *= left;
}

matrix::matrix_line operator/(const matrix::matrix_line& left, double right) {
	matrix::matrix_line res{left};
	return res /= right;
}

void matrix::swap_columns(int i, int j) {
	for (int k = 0; k < n; k++) {
		std::swap((*this)[k][i], (*this)[k][j]);
	}
}

void matrix::swap_rows(int i, int j) {
	std::swap((*this)[i], (*this)[j]);
}

matrix operator+(const matrix& left, const matrix& right) {
	matrix res{left};
	return res += right;
}

matrix operator-(const matrix& left, const matrix& right) {
	matrix res{left};
	return res -= right;
}

matrix operator*(const matrix& left, const matrix& right) {
	matrix res{left};
	return res *= right;
}

matrix operator*(const matrix& left, double right) {
	matrix res{left};
	return res *= right;
} 

matrix operator*(double left, const matrix& right) {
	matrix res{right};
	return res *= left;
}

vector operator*(const matrix& left, const vector& right) {
	vector res{left.height()};

	for (int i = 0; i < left.height(); i++)
		for (int j = 0; j < left.width(); j++)
			res[i] += left.lines[i][j] * right[j];

	return res;
}

matrix& matrix::operator+=(const matrix& m) {
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			lines[i][j] += m.lines[i][j];

	return *this;
}

matrix& matrix::operator-=(const matrix& m) {
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			lines[i][j] -= m.lines[i][j];

	return *this;
}

matrix& matrix::operator*=(const matrix& m) {
	matrix tmp = {this->n}; 
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
				tmp[i][j] += lines[i][k] * m.lines[k][j];

	*this = std::move(tmp);
	return *this;
}

matrix& matrix::operator*=(double d) {
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			lines[i][j] *= d;

	return *this;
}

matrix matrix::inverse() { 
	matrix tmp = *this, res = {n};
	for (int i = 0; i < n; i++)
		res[i][i] = 1;

	for (int i = 0; i < n; i++) {
		int row_max = -1;
		double max = -1;
		for (int j = i; j < n; j++) {
			if (std::abs(lines[j][i]) > max) {
				max = std::abs(lines[j][i]);
				row_max = j;
			}
		}

		res.swap_rows(i, row_max);
		tmp.swap_rows(i, row_max);

		for (int j = i + 1; j < n; j++) {
			double k_j = - tmp[j][i] / tmp[i][i];
			tmp[j] += tmp[i] * k_j;
			res[j] += res[i] * k_j;
		}

		res[i] /= tmp[i][i];
		tmp[i] /= tmp[i][i];
	}

	for (int i = n - 1; i >= 0; i--) {
		for (int j = i - 1; j >= 0; j--) {
			res[j] += res[i] * -tmp[j][i];
			tmp[j] += tmp[i] * -tmp[j][i];
		}
	}

	return res;
}
