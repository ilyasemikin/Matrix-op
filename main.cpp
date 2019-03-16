#include <iostream>
#include <vector>

#include "mx_base.h"



int main() {
	std::vector<long double> mx_items(9, 0);
	
	mx_items = { 0, 1, 2, 3, 4, 5, 6, 7, 9 };;

	Matrix<long double> mx(mx_items, 3, 4, 0.001);
	std::cout << mx;
	mx.toIdentity().print();
	char c;
	std::cin >> c; 
	return 0;
}