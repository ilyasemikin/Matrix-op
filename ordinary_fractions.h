#pragma once

#include <cstdint>

class OrdFract {
private:
	int64_t numerator;
	uint64_t denominator;

public:
	OrdFract(int64_t _numerator = 0, uint64_t _denominator = 1);

	inline bool isExist() const { return (denominator == 0) ? false : true; };
	void reduce();

	OrdFract operator = (const OrdFract &rhs);
	OrdFract operator - ();
	OrdFract operator + ();
	OrdFract operator += (const OrdFract &rhs);
	OrdFract operator -= (const OrdFract &rhs);
	OrdFract operator *= (const OrdFract &rhs);
	OrdFract operator /= (const OrdFract &rhs);

	friend OrdFract operator + (const OrdFract &lhs, const OrdFract &rhs);
	friend OrdFract operator - (const OrdFract &lhs, const OrdFract &rhs);
	friend OrdFract operator * (const OrdFract &lhs, const OrdFract &rhs);
	friend OrdFract operator / (const OrdFract &lhs, const OrdFract &rhs);

	friend std::ostream &operator << (std::ostream &lhs, const OrdFract &rhs);

	friend bool operator == (const OrdFract &lhs, const OrdFract &rhs);
	friend bool operator != (const OrdFract &lhs, const OrdFract &rhs);
	friend bool operator > (const OrdFract &lhs, const OrdFract &rhs);
	friend bool operator >= (const OrdFract &lhs, const OrdFract &rhs);
	friend bool operator < (const OrdFract &lhs, const OrdFract &rhs);
	friend bool operator <= (const OrdFract &lhs, const OrdFract &rhs);
};