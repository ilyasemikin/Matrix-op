#pragma once

#include <vector>

namespace SMath {
	template <typename T>
	T abs(T x);
	template <typename T>
	size_t max(std::vector<T> items);
}

template <typename T>
T SMath::abs(T x) {
	return (x < 0) ? -x : x;
}

template <typename T>
size_t SMath::max(std::vector<T> items) {
	if (items.size() == 0)
		throw(std::exception("Vector is empty"));
	if (items.size() == 1)
		return 0;
	size_t max = 0;
	for (size_t i = 1; i < items.size(); i++)
		if (items[i] > items[max])
			max = i;
	return max;
}