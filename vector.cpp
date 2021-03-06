#include <iostream>
#include <cmath>

#include "vector.h"
#include <vector>

vector::vector(int n) : n(n) {
	data.resize(n);
}

vector::vector(const vector& v) {
	n = v.n;
	data = v.data;
}

vector::vector(vector&& v) {
	n = v.n;
	v.n = 0;
	data = std::move(v.data);
}

vector& vector::operator=(const vector& v) {
	n = v.n;
	data = v.data;

	return *this;
}

vector& vector::operator=(vector&& v) {
	n = v.n;
	v.n = 0;
	data = std::move(v.data);

	return *this;
}

vector::~vector() {
}

double& vector::operator[](int i) {
	return data[i];
}

const double& vector::operator[](int i) const {
	return data[i];
}

void vector::dump() {
	std::cout << "vector [" << this << "]:" << std::endl;
	std::cout << "\tn = " << n << std::endl;
	std::cout << "\tdata = [";
	if (data.empty()) {
		std::cout << "nil]" << std::endl;
	} else {
		for (int i = 0; i < n - 1; i++)
			std::cout << data[i] << ", ";
		std::cout << data[n - 1] << "]" << std::endl;
	}
}

double vector::norm1() {
	double max = 0;
	for (int i = 0; i < n; i++)
		if (max < std::abs(data[i]))
			max = std::abs(data[i]);

	return max;
}

double vector::norm2() {
	double sum = 0;
	for (int i = 0; i < n; i++)
		sum += std::abs(data[i]);

	return sum;
}

double vector::norm3() {
	double sum = 0;
	for (int i = 0; i < n; i++)
		sum += data[i] * data[i];

	return std::sqrt(sum);
}

void vector::swap(int i, int j) {
	std::swap((*this)[i], (*this)[j]);
}

vector operator+(const vector& left, const vector& right) {
	vector res{left};
	return res += right;
}

vector operator-(const vector& left, const vector& right) {
	vector res{left};
	return res -= right;
}

vector operator*(const vector& left, double right) {
	vector res{left};
	return res *= right;
}

vector operator*(double left, const vector& right) {
	vector res{right};
	return res *= left;
}

vector operator/(const vector& left, double right) {
	vector res{left};
	return res /= right;
}

vector& vector::operator+=(const vector& oth) {
	for (int i = 0; i < n; i++)
		this->data[i] += oth.data[i];

	return *this;
}

vector& vector::operator-=(const vector& oth) {
	for (int i = 0; i < n; i++)
		this->data[i] -= oth.data[i];

	return *this;
}

vector& vector::operator*=(double d) {
	for (int i = 0; i < n; i++)
		this->data[i] *= d;

	return *this;
}

vector& vector::operator/=(double d) {
	for (int i = 0; i < n; i++)
		this->data[i] /= d;

	return *this;
}
