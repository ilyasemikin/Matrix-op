#include <iostream>
#include <vector>

#include "mx_base.h"
#include "ordinary_fractions.h"

int main(int argc, char **argv) {
	std::vector<OrdFract> mx_items(9, 0);

	OrdFract of1(-3, 1), of2(3, 1);

	std::cout << of1 / of2 << std::endl;

	mx_items = { 1, 2, 3, 4, 5, 6, 7, 8, 8 };;

	Matrix<OrdFract> mx(mx_items, 3, 3, 0);
	std::cout << mx;
	std::cout << mx.inverse() << std::endl;

	char c;
	std::cin >> c; 
	
	return 0;
}