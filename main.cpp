#include <iostream>
#include <vector>

#include "mx_base.h"
#include "ordinary_fractions.h"

int main(int argc, char **argv) {
	std::vector<long double> mx_items(9, 0),
		mx1_items(9, 0),
		mx2_items(9, 0);
	
	mx_items = { 1, 2, 3, 4, 5, 6, 7, 8, 8 };
	mx1_items = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	mx2_items = { 1, 1, 1, 2, 2, 2, 3, 3, 3 };

	Matrix<long double> mx(mx_items, 3, 3, 0.001),
		mx1(mx1_items, 3, 3, 0.001),
		mx2(mx2_items, 3, 3, 0.001);
	std::cout << mx;

	std::cout << mx1.rank() << std::endl;
	std::cout << mx2.rank() << std::endl;

	std::cout << mx.inverse();

	std::cout << mx.inverseThread() << std::endl;

	char c;
	std::cin >> c; 
	
	return 0;
}