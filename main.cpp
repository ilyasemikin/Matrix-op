#include <iostream>
#include <vector>

#include "mx.h"

int main() {
	std::vector<long double> mx_items(81, 0);
	
	//for (size_t i = 0; i < 10000; i++)
	//	mx_items[i * 100 + i] = 1;
	Matrix<long double> mx(mx_items, 9, 9, 0.001);
	mx.print();
	std::cout << '\n' << mx.determinant() << std::endl;
	char c;
	std::cin >> c; 
	return 0;
}