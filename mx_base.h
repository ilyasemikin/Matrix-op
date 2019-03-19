#pragma once

#include <iostream>
#include <vector>
#include <utility>
#include <iomanip>
#include <sstream>
#include <thread>
#include <algorithm>

#include "smath.h"

template <typename ItemsType>
class Matrix {
private:
	std::vector<ItemsType> items;
	std::size_t m, n;

	std::size_t ij_to_pos(std::size_t i, std::size_t j) const;
	std::pair<std::size_t, std::size_t> pos_to_ij(std::size_t pos) const;

	// epsilon - �������� ���������� ����
	ItemsType epsilon;

	ItemsType _det();
	ItemsType _det2();
	ItemsType _det3();
	// _detn(): ���������� ������������ ������� ������� n ����� �������� ������������ � ������������������ ����
	ItemsType _detn();
	// _crossThread(): ���������� ����������� ������ �������. ������������ ��������������� ������� crossThread()
	void _crossThread(std::size_t posStart, std::size_t posEnd, std::vector<ItemsType> &results);
public:
	Matrix() { m = n = epsilon = 0; };
	Matrix(ItemsType *_items, std::size_t _lenItems, std::size_t _m, std::size_t _n);
	Matrix(std::vector<ItemsType> _items, std::size_t _m, std::size_t _n, ItemsType _epsilon = 0);

	void print(std::ostream &stream) const;
	void print() const;

	void setItem(ItemsType item, std::size_t i, std::size_t j);
	inline void setEpsilon(ItemsType _epsilon) { epsilon = _epsilon; };

	inline std::size_t getM() const { return m; };
	inline std::size_t getN() const { return n; };
	ItemsType getItem(std::size_t i, std::size_t j) const;
	ItemsType getItem(std::size_t pos) const;
	inline ItemsType getEpsilon() const { return epsilon; };

	inline bool isSquare() const { return (m == n && m != 0) ? true : false; };
	inline bool isExist() const { return (m == 0 || n == 0 || items.empty()) ? false : true; }
	// isUpTriangular(): ��������, �������� �� ������� �����������������
	bool isUpTriangular() const;
	// isUpTriangular(): ��������, �������� �� ������� ����������������
	bool isLowTriangular() const;
	// isZeroLine(i): ��������, ������� �� i-������ �� �����
	bool isZeroLine(std::size_t i) const;
	// isZeroColumn(j): ��������, ������� �� j-������� �� �����
	bool isZeroColumn(std::size_t j) const;

	void insertLine(std::size_t i, std::vector<ItemsType> line);
	void insertColumn(std::size_t j, std::vector<ItemsType> column);

	void addLine(std::size_t i, std::vector<ItemsType> line);
	void addColumn(std::size_t j, std::vector<ItemsType> column);

	// getLineV(i): ���������� i ������ � ���� vector<��� �������� �������>
	std::vector<ItemsType> getLineV(std::size_t i);
	// getColumn(j): ���������� j ������� � ���� vector<��� �������� �������>
	std::vector<ItemsType> getColumnV(std::size_t j);
	// getLine(i): ���������� i ������ � ���� Matrixr<��� �������� �������>
	Matrix<ItemsType> getLine(std::size_t i);
	// getColumn(j): ���������� j ������� � ���� Matrix<��� �������� �������>
	Matrix<ItemsType> getColumn(std::size_t j);

	void swapLines(std::size_t i1, std::size_t i2);
	void swapColumns(std::size_t j1, std::size_t j2);

	// minor(...): ���������� ����� �������, �� �������� �������� �������
	Matrix<ItemsType> minor(std::size_t pos) const;
	Matrix<ItemsType> minor(std::size_t _i, std::size_t _j) const;
	// transposed(): ������� ����������������� �������, �� �������� �������� �������
	Matrix<ItemsType> transposed();
	// cross(): ������� �������� �������, �� �������� �������� �������
	Matrix<ItemsType> cross();
	// crossThread(): ������� �������� �������, �� �������� �������� �������. ���������� � �������������� ���������������
	Matrix<ItemsType> crossThread();
	// inverse(): ������� �������� �������, �� �������� �������� �������
	Matrix<ItemsType> inverse();
	// inverseThread(): ������� �������� �������, �� �������� �������� �������. ���������� � �������������� ���������������
	Matrix<ItemsType> inverseThread();

	// toRowEchelonForm(): �������� ������� � ������ ������������ ����
	Matrix<ItemsType> toRowEchelonForm();
	// toIdentity(): �������� ������� � ���������
	Matrix<ItemsType> toIdentity();

	// determinant(): ���������� ������������ �������
	ItemsType determinant();
	// rank(): ����������� ����� �������
	std::size_t rank();

	Matrix<ItemsType> operator = (const Matrix<ItemsType> &mx);
	Matrix<ItemsType> operator + () const;
	Matrix<ItemsType> operator - () const;

	template <typename T>
	friend Matrix<T> operator + (const Matrix<T> &mx1, const Matrix<T> &mx2);
	template <typename T, typename NumType>
	friend Matrix<T> operator + (const Matrix<T> &mx, const NumType &num);
	template <typename T, typename NumType>
	friend Matrix<T> operator + (const NumType &num, const Matrix<T> &mx);
	template <typename T>
	friend Matrix<T> operator - (const Matrix<T> &mx1, const Matrix<T> &mx2);
	template <typename T, typename NumType>
	friend Matrix<T>  operator - (const Matrix<T> &mx, const NumType &nun);
	template <typename T, typename NumType>
	friend Matrix<T> operator - (const NumType &num, const Matrix<T> &mx);
	template <typename T>
	friend Matrix<T> operator * (const Matrix<T> &mx1, const Matrix<T> &mx2);
	template <typename T, typename NumType>
	friend Matrix<T> operator * (const Matrix<T> &mx, const NumType &num);
	template <typename T, typename NumType>
	friend Matrix<T> operator * (const NumType &num, const Matrix<T> &mx);

	template <typename T>
	friend bool operator == (const Matrix<T> &mx1, const Matrix<T> &mx2);

	template <typename T>
	friend std::ostream & operator << (std::ostream &stream, const Matrix<ItemsType> &mx);
};

template <typename ItemsType>
std::size_t Matrix<ItemsType>::ij_to_pos(std::size_t i, std::size_t j) const {
	return i * n + j;
}

template <typename ItemsType>
std::pair<std::size_t, std::size_t> Matrix<ItemsType>::pos_to_ij(std::size_t pos) const {
	std::size_t i, j;
	i = pos / n;
	j = pos % n;
	return std::pair<std::size_t, std::size_t>(i, j);
}

// _det(): ���������� ������� ����� ���������� ������ ������. �� ������ ������ �� ������������. ��������� ���� �� ��� �������, ��� ����� ������� (����� ����������)
template <typename ItemsType>
ItemsType Matrix<ItemsType>::_det() {
	ItemsType res = 0;
	if (m == 1)
		return items[0];
	else if (m == 2)
		return _det2();
	else if (m == 3)
		return _det3();
	else {
		ItemsType term;
		std::size_t i;
		for (i = 0; i < m; i++) {
			term = (((i % 2) == 0) ? 1 : -1) * items[i] * minor(i)._det();
			res += term;
		}
	}
	return res;
}

template <typename ItemsType>
ItemsType Matrix<ItemsType>::_det2() {
	return items[0] * items[3] - items[1] * items[2];
}

template <typename ItemsType>
ItemsType Matrix<ItemsType>::_det3() {
	return items[0] * items[4] * items[8] + items[1] * items[5] * items[6] + items[2] * items[3] * items[7] - items[2] * items[4] * items[6] - items[1] * items[3] * items[8] - items[0] * items[5] * items[7];
}

template <typename ItemsType>
ItemsType Matrix<ItemsType>::_detn() {
	ItemsType det = 1,
		offsetCoef = 1;			// ����������, ������� ��������� ��� ������� ������ �� �����, ������� ������������ ��� ��������� � ������� i, i ������� (�.� �����, ������� ���������)
	bool sign = false;			// ���� ������������, �������� �� �� ����, ��� � �������� ��������, � ������� (i == j) ����� �������� ������
	Matrix<ItemsType> mx = *this;

	for (std::size_t j = 0; j < mx.m; j++) {
		for (std::size_t i = j; i < m; i++) {
			if (i == j) {
				// ��������� ������� ���������� ����, ����� � ������� �� ������� i, i ��������� �������
				if (SMath::abs(mx.items[mx.ij_to_pos(i, j)] - 1) > epsilon)
					if (SMath::abs(mx.items[mx.ij_to_pos(i, j)]) <= epsilon) {
						std::vector<ItemsType> curColumn = mx.getColumnV(j);
						curColumn.erase(curColumn.begin(), curColumn.begin() + i);
						std::size_t iMaxItem = SMath::max(curColumn) + i;
						if (SMath::abs(mx.items[mx.ij_to_pos(iMaxItem, j)]) <= epsilon)
							return 0;

						mx.swapLines(i, iMaxItem);
						std::cout << "Swap lines: " << i << ' ' << iMaxItem << std::endl;
						sign = !sign;
						if (SMath::abs(mx.items[mx.ij_to_pos(i, j)] - 1) <= epsilon)
							continue;
					}

				ItemsType coef = mx.items[mx.ij_to_pos(i, j)];
				offsetCoef *= coef;
				for (std::size_t k = j; k < mx.m; k++)
					mx.items[mx.ij_to_pos(i, k)] /= coef;

				continue;
			}
			// ������ ������� � �����������������
			ItemsType coef = mx.items[mx.ij_to_pos(i, j)];
			for (std::size_t k = j; k < mx.m; k++)
				mx.items[mx.ij_to_pos(i, k)] -= coef * mx.items[mx.ij_to_pos(j, k)];
		}
	}

	// ������� ������������ ����������������� �������
	for (std::size_t i = 0; i < m; i++)
		det *= mx.items[ij_to_pos(i, i)];

	return (sign ? -1 : 1) * offsetCoef * det;
}

template <typename ItemsType>
void Matrix<ItemsType>::_crossThread(std::size_t posStart, std::size_t posEnd, std::vector<ItemsType> &results) {
	for (size_t i = posStart; i < posEnd; i++)
		results[i] = ((i % 2 == 0) ? 1 : -1) * (this->minor(i)).determinant();
}

template <typename ItemsType>
Matrix<ItemsType>::Matrix(ItemsType *_items, std::size_t _lenItems, std::size_t _m, std::size_t _n) {
	items.resize(_m * _n);
	for (std::size_t i = 0; i < _m * _n; i++)
		items[i] = (i < _lenItems) ? _items[i] : 0;
	m = _m, n = _n;
	epsilon = 0;
}

template <typename ItemsType>
Matrix<ItemsType>::Matrix(std::vector<ItemsType> _items, std::size_t _m, std::size_t _n, ItemsType _epsilon) {
	items.resize(_m * _n);
	for (std::size_t i = 0; i < _m * _n; i++)
		items[i] = (i < _items.size()) ? _items[i] : 0;
	m = _m, n = _n;
	epsilon = _epsilon;
}

template <typename ItemsType>
ItemsType Matrix<ItemsType>::getItem(std::size_t i, std::size_t j) const {
	if (i >= m || j >= n)
		throw std::exception("Item not defined");
	return items[ij_to_pos(i, j)];
}

template <typename ItemsType>
ItemsType Matrix<ItemsType>::getItem(std::size_t pos) const {
	if (pos >= m * n)
		throw std::exception();
	return items[pos];
}

template <typename ItemsType>
void Matrix<ItemsType>::print(std::ostream &stream) const {
	if (!isExist())
		throw(std::exception("Matrix not defined"));
	std::size_t sizePlaceItem = 1;			// �����, ��������� ��� ����� � �������
	for (auto &x : items) {
		std::stringstream itemStream;
		itemStream << std::setprecision(2) << x;
		std::size_t itemLength = itemStream.str().size();
		if (sizePlaceItem < itemLength)
			sizePlaceItem = itemLength;
	}
	for (std::size_t i = 0; i < m; i++) {
		stream << "| ";
		for (std::size_t j = 0; j < n; j++)
			stream << std::setw(sizePlaceItem) << std::setprecision(2) << items[ij_to_pos(i, j)] << ' ';
		stream << '|' << std::endl;
	}
}

template <typename ItemsType>
void Matrix<ItemsType>::print() const {
	print(std::cout);
}

template <typename ItemsType>
void Matrix<ItemsType>::setItem(ItemsType item, std::size_t i, std::size_t j) {
	if (!isExist() || i >= m || j >= n)
		throw(std::exception("Position item not defined"));
	items[i, j] = item;
}

template <typename ItemsType>
bool Matrix<ItemsType>::isUpTriangular() const {
	for (std::size_t i = 0; i < m; i++)
		for (std::size_t j = i + 1; j < n; j++)
			if (!(SMath::abs(items[ij_to_pos(i, j)]) <= epsilon))
				return false;
	return true;
}

template <typename ItemsType>
bool Matrix<ItemsType>::isLowTriangular() const {
	for (std::size_t i = m - 1; i > 0; i--)
		for (std::size_t j = 0; j < m - i - 1; j++)
			if (!(items[ij_to_pos(i, j)] <= 0))
				return false;
	return true;
}

template <typename ItemsType>
bool Matrix<ItemsType>::isZeroLine(std::size_t i) const {
	for (std::size_t k = 0; k < n; k++)
		if (SMath::abs(items[ij_to_pos(i, k)]) > epsilon)
			return false;
	return true;
}

template <typename ItemsType>
bool Matrix<ItemsType>::isZeroColumn(std::size_t j) const {
	for (std::size_t k = 0; k < m; k++)
		if (SMath::abs(items[ij_to_pos(k, j)]) > epsilon)
			return false;
	return true;
}

template <typename ItemsType>
void Matrix<ItemsType>::insertLine(std::size_t i, std::vector<ItemsType> line) {
	if (i >= m)
		throw std::exception("Line not defined");
	std::size_t k = 0;
	for (std::size_t j = 0; j < n; j++)
		items[ij_to_pos(i, j)] = (k < line.size()) ? line[k++] : 0;
}

template <typename ItemsType>
void Matrix<ItemsType>::insertColumn(std::size_t j, std::vector<ItemsType> column) {
	if (j >= n)
		throw std::exception("Column not defined");
	std::size_t k = 0;
	for (std::size_t i = 0; i < m; i++)
		items[ij_to_pos(i, j)] = (k < column.size()) ? column[k++] : 0;
}

template <typename ItemsType>
void Matrix<ItemsType>::addLine(std::size_t i, std::vector<ItemsType> line) {
	m = (i <= m) ? m + 1 : i + 1;
	items.resize(m * n, 0);
	std::size_t k = 0;
	for (std::size_t j = 0; j < n; j++)
		items.insert(items.begin() + ij_to_pos(i, j), (k < line.size()) ? line[k++] : 0);
}

template <typename ItemsType>
void Matrix<ItemsType>::addColumn(std::size_t j, std::vector<ItemsType> column) {
	n = (j <= n) ? n + 1 : j + 1;
	items.resize(m * n, 0);
	std::size_t k = 0;
	for (std::size_t i = 0; i < m; i++)
		items.insert(items.begin() + ij_to_pos(i, j), (k < column.size()) ? column[k++] : 0);
}

template <typename ItemsType>
std::vector<ItemsType> Matrix<ItemsType>::getLineV(std::size_t i) {
	if (i >= m)
		throw(std::exception("Line not defined"));
	std::vector<ItemsType> res(n);
	for (std::size_t j = 0; j < n; j++)
		res[j] = items[ij_to_pos(i, j)];
	return res;
}

template <typename ItemsType>
std::vector<ItemsType> Matrix<ItemsType>::getColumnV(std::size_t j) {
	if (j >= n)
		throw(std::exception("Column not defined"));
	std::vector<ItemsType> res(m);
	for (std::size_t i = 0; i < m; i++)
		res[i] = items[ij_to_pos(i, j)];
	return res;
}

template <typename ItemsType>
Matrix<ItemsType> Matrix<ItemsType>::getLine(std::size_t i) {
	return Matrix<ItemsType>(getLineV(i), 1, n);
}

template <typename ItemsType>
Matrix<ItemsType> Matrix<ItemsType>::getColumn(std::size_t j) {
	return Matrix<ItemsType>(getColumnV(j), m, 1);
}

template <typename ItemsType>
void Matrix<ItemsType>::swapLines(std::size_t i1, std::size_t i2) {
	if (i1 >= m || i2 >= m)
		throw(std::exception("Lines not defined"));
	for (std::size_t j = 0; j < n; j++) {
		ItemsType temp;
		temp = items[ij_to_pos(i1, j)];
		items[ij_to_pos(i1, j)] = items[ij_to_pos(i2, j)];
		items[ij_to_pos(i2, j)] = temp;
	}
}

template <typename ItemsType>
void Matrix<ItemsType>::swapColumns(std::size_t j1, std::size_t j2) {
	if (j1 >= n || j2 >= n)
		throw(std::exception("Lines not defined"));
	for (std::size_t i = 0; i < m; i++) {
		ItemsType temp;
		temp = items[ij_to_pos(i, j1)];
		items[ij_to_pos(i, j1)] = items[ij_to_pos(i, j2)];
		items[ij_to_pos(i, j2)] = temp;
	}
}

template <typename ItemsType>
Matrix<ItemsType> Matrix<ItemsType>::minor(std::size_t pos) const {
	std::pair<std::size_t, std::size_t> ij = pos_to_ij(pos);
	return minor(ij.first, ij.second);
}

template <typename ItemsType>
Matrix<ItemsType> Matrix<ItemsType>::minor(std::size_t _i, std::size_t _j) const {
	if (_i >= m || _j >= n || m == 1 || n == 1)
		throw std::exception("Minor not defined");
	std::vector<ItemsType> minorItems((m - 1) * (n - 1));
	std::size_t k = 0;
	for (std::size_t i = 0; i < m; i++) {
		if (i == _i)
			continue;
		for (std::size_t j = 0; j < n; j++) {
			if (j == _j)
				continue;
			minorItems[k++] = items[ij_to_pos(i, j)];
		}
	}
	return Matrix<ItemsType>(minorItems, m - 1, n - 1);
}

template <typename ItemsType>
Matrix<ItemsType> Matrix<ItemsType>::transposed() {
	std::vector<ItemsType> transItems(m * n);
	for (std::size_t i = 0; i < m; i++)
		for (std::size_t j = 0; j < n; j++)
			transItems[ij_to_pos(i, j)] = items[ij_to_pos(j, i)];
	return Matrix<ItemsType>(transItems, m, n);
}

template <typename ItemsType>
Matrix<ItemsType> Matrix<ItemsType>::cross() {
	std::vector<ItemsType> crossItems(m * n);
	for (std::size_t i = 0; i < m; i++)
		for (std::size_t j = 0; j < n; j++)
			crossItems[ij_to_pos(i, j)] = ((((i + j) % 2) == 0) ? 1 : -1) * minor(i, j).determinant();
	return Matrix<ItemsType>(crossItems, m, n).transposed();
}

template <typename ItemsType>
Matrix<ItemsType> Matrix<ItemsType>::crossThread() {
	std::vector<ItemsType> results(m * n);		// �������������� ������
	std::size_t const length = std::distance(items.begin(), items.end());			// ���������� ��������� �� ���������
	std::size_t const minPerThread = 2;												// ����������� ���������� ��������� �� ��������� �� ����� ������
	std::size_t const maxThreads = (length * minPerThread - 1) / minPerThread;		// ���������� �������� ���������� �������
	std::size_t const hardwareThreads = std::thread::hardware_concurrency();		// ���������� �������, �������������� ��������
	std::size_t const numThreads = std::min((hardwareThreads != 0) ? hardwareThreads : 2, maxThreads);	// ���������� �������, ������� ����� ������������� �������������
	std::size_t const blockSize = length / numThreads;								// ���������� ���������, �������������� ����� �������
	std::vector<std::thread> threads(numThreads - 1);								// ������� �������
	std::size_t itemStart = 0;
	for (std::size_t i = 0; i < (numThreads - 1); i++) {
		std::size_t itemEnd = itemStart + blockSize;
		threads[i] = std::thread(&Matrix<ItemsType>::_crossThread, this, itemStart, itemEnd, std::ref(results));		// �������� �������	
		std::size_t itemStart = itemEnd;
	}
	_crossThread(itemStart, items.size(), std::ref(results));
	std::for_each(threads.begin(), threads.end(), std::mem_fn(&std::thread::join));		// ����������� � ��������� �������
	return Matrix<ItemsType>(results, m, n).transposed();
}

template <typename ItemsType>
Matrix<ItemsType> Matrix<ItemsType>::inverse() {
	ItemsType det = determinant();
	if (SMath::abs(det) <= epsilon)
		throw(std::exception("Inverse matrix is not defined"));
	std::vector<ItemsType> inverseItems(m * n);
	for (std::size_t i = 0; i < m; i++)
		for (std::size_t j = 0; j < n; j++)
			inverseItems[ij_to_pos(i, j)] = ((((i + j) % 2) == 0) ? 1 : -1) * minor(i, j).determinant() / det;
	return Matrix<ItemsType>(inverseItems, m, n).transposed();
}

template <typename ItemsType>
Matrix<ItemsType> Matrix<ItemsType>::inverseThread() {
	Matrix<ItemsType> res = crossThread();
	ItemsType det = determinant();
	for (auto &x : res.items)
		x /= det;
	return res;
}

template <typename ItemsType>
Matrix<ItemsType> Matrix<ItemsType>::toRowEchelonForm() {
	if (!isExist())
		throw(std::exception("Matrix isn't exist"));
	if (m > n)
		throw;
	if (isLowTriangular())
		return *this;
	Matrix<ItemsType> res = *this;
	for (std::size_t j = 0; j < res.m; j++)
		for (std::size_t i = 0; i < res.m; i++) {
			if (i == j) {
				if (SMath::abs(res.items[ij_to_pos(i, j)]) <= res.epsilon) {
					std::vector<ItemsType> column = res.getColumnV(j);
					column.erase(column.begin(), column.begin() + j);
					std::size_t iMaxItem = SMath::max(column) + j;
					if (SMath::abs(items[ij_to_pos(iMaxItem, j)]) < res.epsilon)
						break;
					swapLines(iMaxItem, i);
				}
				continue;
			}

			if (i > j) {
				ItemsType coef0, coef1;
				coef0 = res.items[ij_to_pos(j, j)];
				coef1 = res.items[ij_to_pos(i, j)];
				for (std::size_t k = j; k < res.n; k++)
					res.items[ij_to_pos(i, k)] = coef0 * res.items[ij_to_pos(i, k)] - coef1 * res.items[ij_to_pos(j, k)];
			}
		}
	return res;
}

template <typename ItemsType>
Matrix<ItemsType> Matrix<ItemsType>::toIdentity() {
	if (!isExist())
		throw(std::exception("Matrix isn't exist"));
	if (m > n)
		throw;
	Matrix<ItemsType> res = *this;
	for (std::size_t j = 0; j < res.m; j++)
		for (std::size_t i = j; i < res.m; i++) {
			if (i == j) {				// �������� �������� i, i � �������
				if (SMath::abs(res.items[ij_to_pos(i, j)]) <= res.epsilon) {
					std::vector<ItemsType> curColumn = res.getColumnV(j);
					curColumn.erase(curColumn.begin(), curColumn.begin() + i);
					std::size_t iMaxItem = SMath::max(curColumn) + i;
					if (SMath::abs(res.items[iMaxItem]) < res.epsilon)
						break;
					res.swapLines(i, iMaxItem);
				}

				if (SMath::abs(res.items[ij_to_pos(i, j)] - 1) > res.epsilon) {
					ItemsType coef = res.items[ij_to_pos(i, j)];
					for (std::size_t k = i; k < res.n; k++)
						res.items[ij_to_pos(i, k)] /= coef;
				}
				continue;
			}

			// �������� ������� � ������ ������������ ����
			ItemsType coef0, coef1;
			coef0 = res.items[ij_to_pos(j, j)];
			coef1 = res.items[ij_to_pos(i, j)];
			for (std::size_t k = j; k < res.n; k++)
				res.items[ij_to_pos(i, k)] = coef0 * res.items[ij_to_pos(i, k)] - coef1 * res.items[ij_to_pos(j, k)];
		}

	// �������� ������� � ���������
	for (std::size_t j = res.m - 1; j > 0; j--)
		for (std::size_t i = 0; i != j; i++) {
			ItemsType coef = res.items[ij_to_pos(i, j)];
			for (std::size_t k = j; k < res.n; k++)
				res.items[ij_to_pos(i, k)] = res.items[ij_to_pos(i, k)] - coef * res.items[ij_to_pos(j, k)];
		}
	return res;
}

template <typename ItemsType>
ItemsType Matrix<ItemsType>::determinant() {
	if (!isExist() || !isSquare())
		throw std::exception("Determinant not defined");
	if (n == 1)
		return items[0];
	else if (n == 2)
		return _det2();
	else if (n == 3)
		return _det3();
	return _detn();
}

template <typename ItemsType>
std::size_t Matrix<ItemsType>::rank() {
	Matrix<ItemsType> rowEchelon = toRowEchelonForm();
	std::size_t res = n;
	for (size_t i = 0; i < n; i++)
		if (rowEchelon.isZeroLine(i))
			res--;
	return res;
}

template <typename ItemsType>
Matrix<ItemsType> Matrix<ItemsType>::operator = (const Matrix<ItemsType> &mx) {
	items = mx.items;
	m = mx.m;
	n = mx.n;
	epsilon = mx.epsilon;
	return Matrix<ItemsType>(mx.items, mx.m, mx.n, mx.epsilon);
}

template <typename ItemsType>
Matrix<ItemsType> Matrix<ItemsType>::operator + () const {
	return this;
}

template <typename ItemsType>
Matrix<ItemsType> Matrix<ItemsType>::operator - () const {
	std::vector<ItemsType> resItems = items;
	for (std::size_t i = 0; i < m * n; i++) {
		if (resItems[i] <= epsilon)
			continue;
		resItems[i] = -resItems[i];
	}
	return Matrix<ItemsType>(resItems, m, n);
}

template <typename ItemsType>
Matrix<ItemsType> operator + (const Matrix<ItemsType> &mx1, const Matrix<ItemsType> &mx2) {
	if (!mx1.isExist() || mx1.m != mx2.m || mx1.n != mx2.n)
		throw(std::exception("Amount not defined"));
	std::vector<ItemsType> resItems(mx1.m * mx1.n);
	for (std::size_t i = 0; i < mx1.m * mx1.n; i++)
		resItems[i] = mx1.items[i] + mx2.items[i];
	return Matrix<ItemsType>(resItems, mx1.m, mx1.n);
}

template <typename ItemsType, typename NumType>
Matrix<ItemsType> operator + (const Matrix<ItemsType> &mx, const NumType &num) {
	if (!mx.isSquare())
		throw("Amount not defined");
	std::vector<ItemsType> resItems = mx.items;
	for (std::size_t i = 0; i < mx.m; i++)
		resItems[mx.ij_to_pos(i, i)] = num + mx.items[mx.ij_to_pos(i, i)];
	return Matrix<ItemsType>(resItems, mx.n, mx.n);
}

template <typename ItemsType, typename NumType>
Matrix<ItemsType> operator + (const NumType &num, const Matrix<ItemsType> &mx) {
	return mx + num;
}

template <typename ItemsType>
Matrix<ItemsType> operator - (const Matrix<ItemsType> &mx1, const Matrix<ItemsType> &mx2) {
	if (!mx1.isExist() || mx1.m != mx2.m || mx1.n != mx2.n)
		throw(std::exception("Difference not defined"));
	return mx1 + (-mx2);
}

template <typename ItemsType, typename NumType>
Matrix<ItemsType> operator - (const Matrix<ItemsType> &mx, const NumType &num) {
	return mx + (-num);
}

template <typename ItemsType, typename NumType>
Matrix<ItemsType> operator - (const NumType &num, const Matrix<ItemsType> &mx) {
	return num + (-mx);
}

template <typename ItemsType>
Matrix<ItemsType> operator * (const Matrix<ItemsType> &mx1, const Matrix<ItemsType> &mx2) {
	if (!mx1.isExist() || mx1.n != mx2.m)
		throw(std::exception("Multiply not defined"));
	std::vector<ItemsType> res(mx1.m * mx2.n);
	for (std::size_t i = 0; i < mx1.m; i++)
		for (std::size_t j = 0; j < mx2.n; j++) {
			res[i * mx2.n + j] = 0;
			for (std::size_t k = 0; k < mx1.n; k++)
				res[i * mx2.n + j] += mx1.items[i * mx1.n + k] * mx2.items[k * mx2.n + j];
		}
	return Matrix<ItemsType>(res, mx1.m, mx2.n);
}

template <typename ItemsType, typename NumType>
Matrix<ItemsType> operator * (const Matrix<ItemsType> &mx, const NumType &num) {
	std::vector<ItemsType> res(mx.m * mx.n);
	for (std::size_t i = 0; i < mx.m * mx.n; i++)
		res[i] = num * mx.items[i];
	return Matrix<ItemsType>(res, mx.m, mx.n);
}

template <typename T, typename NumType>
Matrix<T> operator * (const NumType &num, const Matrix<T> &mx) {
	return mx * num;
}

template <typename ItemsType>
bool operator == (const Matrix<ItemsType> &mx1, const Matrix<ItemsType> &mx2) {
	if (!mx1.isExist() || !mx2.isExist() || mx1.m != mx2.m || mx1.n != mx2.n)
		return false;
	for (std::size_t i = 0; i < mx1.m * mx1.n; i++)
		if (mx1.items[i] != mx2.items[i])
			return false;
	return true;
}

template<typename ItemsType>
std::ostream & operator<<(std::ostream &stream, const Matrix<ItemsType>&mx) {
	mx.print(stream);
	return stream;
}