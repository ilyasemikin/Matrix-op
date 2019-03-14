#include <iostream>
#include <vector>

#include "mx.h"



int main() {
	std::vector<long double> mx_items(9, 0);
	
	mx_items = { 0, 1, 2, 3, 4, 5, 6, 7, 9 };;

	//for (size_t i = 0; i < 10000; i++)
	//	mx_items[i * 100 + i] = 1;
	Matrix<long double> mx(mx_items, 3, 3, 0.001);
	mx.print();
	std::cout << '\n' << mx.determinant() << std::endl;
	char c;
	std::cin >> c; 
	return 0;
}