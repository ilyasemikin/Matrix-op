#include <iostream>
#include <vector>

#include "mx_base.h"
#include "ordinary_fractions.h"

int main(int argc, char **argv) {
	std::vector<long double> mx_items(9, 0),
		mx1_items(9, 0),
		mx2_items(9, 0);
	mx_items = { 2, -1, -1, 4, 5, 3, 4, -2, 11, 8, 3, -2, 4, 11, 3, 1, 1, 1, 1, 5 };
	mx1_items = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	mx2_items = { 1, 1, 1, 2, 2, 2, 3, 3, 3 };

	Matrix<long double> mx(mx_items, 6, 5, 0.001);

	std::cout << mx << std::endl;

	std::cout << mx.toRowEchelonForm() << std::endl;

	std::cout << mx.toIdentity() << std::endl;

	system("pause");
	return 0;
}