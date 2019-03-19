#include <iostream>
#include <vector>

#include "mx_base.h"
#include "ordinary_fractions.h"

int main(int argc, char **argv) {
	std::vector<long double> mx_items(9, 0);

	mx_items = { 1, 2, 3, 4, 5, 6, 7, 8, 8 };

	Matrix<long double> mx(mx_items, 3, 3, 0);

	std::cout << mx;

	std::cout << mx.inverse();

	std::cout << mx.inverseThread() << std::endl;

	char c;
	std::cin >> c; 
	
	return 0;
}