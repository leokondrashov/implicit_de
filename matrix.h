#ifndef MATRIX_H
#define MATRIX_H

#include <vector>

#include "vector.h"

class matrix {
public:
	class matrix_line {
	private:
		int n;
		std::vector<double> data;

	public:
		matrix_line(int n);
		matrix_line(const matrix_line& ml);
		matrix_line(matrix_line&& ml);

		matrix_line& operator=(const matrix_line& ml);
		matrix_line& operator=(matrix_line&& ml);

		~matrix_line();

		double& operator[](int i);
		const double& operator[](int i) const;

		matrix_line& operator+=(const matrix_line& ml);
		matrix_line& operator*=(double d);
		matrix_line& operator/=(double d);

		friend matrix_line operator+(const matrix_line& left, const matrix_line& right);
		friend matrix_line operator*(const matrix_line& left, double right);
		friend matrix_line operator*(double left, const matrix_line& right);
		friend matrix_line operator/(const matrix_line& left, double right);

	};

private:
	int n;
	std::vector<matrix_line> lines;

public:
	matrix(int n);
	matrix(const matrix& m);
	matrix(matrix&& m);

	matrix& operator=(const matrix& m);
	matrix& operator=(matrix&& m);

	~matrix();

	matrix_line& operator[](int i);
	const matrix_line& operator[](int i) const;

	friend matrix operator+(const matrix& left, const matrix& right);
	friend matrix operator-(const matrix& left, const matrix& right);
	friend matrix operator*(const matrix& left, const matrix& right);
	friend matrix operator*(const matrix& left, double right);
	friend matrix operator*(double left, const matrix& right);
	friend vector operator*(const matrix& left, const vector& right);
	
	matrix& operator+=(const matrix& m);
	matrix& operator-=(const matrix& m);
	matrix& operator*=(const matrix& m);
	matrix& operator*=(double d);

	matrix inverse();

	void swap_columns(int i, int j);
	void swap_rows(int i, int j);

	void dump();

	int width() const { return n; }
	int height() const { return n; }
};

#endif // MATRIX_H

