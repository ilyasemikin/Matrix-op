#include <numeric>
#include <iostream>
#include "ordinary_fractions.h"
#include "smath.h"

OrdFract::OrdFract(int64_t _numerator, uint64_t _denominator) {
	numerator = _numerator;
	denominator = _denominator;
	reduce();
}

void OrdFract::reduce() {
	uint64_t divider = std::gcd(SMath::abs(numerator), denominator);
	numerator /= int64_t(divider);
	denominator /= divider;
}

OrdFract OrdFract::operator = (const OrdFract &rhs) {
	numerator = rhs.numerator;
	denominator = rhs.denominator;
	reduce();
	return OrdFract(numerator, denominator);
}

OrdFract OrdFract::operator - () {
	return OrdFract(-numerator, denominator);
}

OrdFract OrdFract::operator + () {
	return *this;
}

OrdFract OrdFract::operator += (const OrdFract &rhs) {
	*this = *this + rhs;
	return *this;
}

OrdFract OrdFract::operator -= (const OrdFract &rhs) {
	*this = *this - rhs;
	return *this;
}

OrdFract OrdFract::operator *= (const OrdFract &rhs) {
	*this = *this * rhs;
	return *this;
}

OrdFract OrdFract::operator /= (const OrdFract &rhs) {
	*this = *this / rhs;
	return *this;
}

OrdFract operator + (const OrdFract &lhs, const OrdFract &rhs) {
	return OrdFract(lhs.numerator * rhs.denominator + rhs.numerator * lhs.denominator, lhs.denominator * rhs.denominator);
}

OrdFract operator - (const OrdFract &lhs, const OrdFract &rhs) {
	return OrdFract(lhs.numerator * rhs.denominator - rhs.numerator * lhs.denominator, lhs.denominator * rhs.denominator);;
}

OrdFract operator * (const OrdFract &lhs, const OrdFract &rhs) {
	return OrdFract(lhs.numerator * rhs.numerator, lhs.denominator * rhs.denominator);
}

OrdFract operator / (const OrdFract &lhs, const OrdFract &rhs) {
	return OrdFract(lhs.numerator * int64_t(rhs.denominator), uint64_t(SMath::abs(rhs.numerator)) * lhs.denominator);
}

std::ostream &operator << (std::ostream &lhs, const OrdFract &rhs) {
	if (!rhs.isExist())
		throw;
	lhs << rhs.numerator;
	if (rhs.denominator != 1)
		lhs << '/' << rhs.denominator;
	return lhs;
}

bool operator == (const OrdFract &lhs, const OrdFract &rhs) {
	return (lhs.numerator * rhs.denominator == lhs.denominator * rhs.numerator);
}

bool operator != (const OrdFract &lhs, const OrdFract &rhs) {
	return (!(lhs == rhs));
}

bool operator > (const OrdFract &lhs, const OrdFract &rhs) {
	return (lhs.numerator * rhs.denominator > lhs.denominator * rhs.numerator);
}

bool operator >= (const OrdFract &lhs, const OrdFract &rhs) {
	return ((lhs > rhs) && (lhs == rhs));
}

bool operator < (const OrdFract &lhs, const OrdFract &rhs) {
	return !(lhs >= rhs);
}

bool operator <= (const OrdFract &lhs, const OrdFract &rhs) {
	return !(lhs > rhs);
}