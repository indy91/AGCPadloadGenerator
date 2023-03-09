// Copyright (c) Martin Schweiger
// Licensed under the MIT License

#pragma once

#include <math.h>

const double PI = 3.14159265358979323846;	///< pi
const double PI05 = 1.57079632679489661923;	///< pi/2
const double PI2 = 6.28318530717958647693;	///< pi*2
const double RAD = PI / 180.0;				///< factor to map degrees to radians
const double DEG = 180.0 / PI;				///< factor to map radians to degrees
const double C0 = 299792458.0;				///< speed of light in vacuum [m/s]
const double TAUA = 499.004783806;			///< light time for 1 AU [s]
const double AU = C0 * TAUA;				///< astronomical unit (mean geocentric distance of the sun) [m]

typedef union {
	double data[3];               ///< array data interface
	struct { double x, y, z; };   ///< named data interface
} VECTOR3;

typedef union {
	double data[9];               ///< array data interface (row-sorted)
	struct { double m11, m12, m13, m21, m22, m23, m31, m32, m33; }; ///< named data interface
} MATRIX3;

inline MATRIX3 _M(double m11, double m12, double m13,
	double m21, double m22, double m23,
	double m31, double m32, double m33)
{
	MATRIX3 mat = { m11,m12,m13,  m21,m22,m23,  m31,m32,m33 };
	return mat;
}

inline MATRIX3 mul(const MATRIX3 &A, const MATRIX3 &B)
{
	MATRIX3 mat = {
		A.m11*B.m11 + A.m12*B.m21 + A.m13*B.m31, A.m11*B.m12 + A.m12*B.m22 + A.m13*B.m32, A.m11*B.m13 + A.m12*B.m23 + A.m13*B.m33,
		A.m21*B.m11 + A.m22*B.m21 + A.m23*B.m31, A.m21*B.m12 + A.m22*B.m22 + A.m23*B.m32, A.m21*B.m13 + A.m22*B.m23 + A.m23*B.m33,
		A.m31*B.m11 + A.m32*B.m21 + A.m33*B.m31, A.m31*B.m12 + A.m32*B.m22 + A.m33*B.m32, A.m31*B.m13 + A.m32*B.m23 + A.m33*B.m33
	};
	return mat;
}

inline VECTOR3 _V(double x, double y, double z)
{
	VECTOR3 vec = { x,y,z }; return vec;
}

inline VECTOR3 mul(const MATRIX3 &A, const VECTOR3 &b)
{
	return _V(
		A.m11*b.x + A.m12*b.y + A.m13*b.z,
		A.m21*b.x + A.m22*b.y + A.m23*b.z,
		A.m31*b.x + A.m32*b.y + A.m33*b.z);
}

inline VECTOR3 tmul(const MATRIX3 &A, const VECTOR3 &b)
{
	return _V(
		A.m11*b.x + A.m21*b.y + A.m31*b.z,
		A.m12*b.x + A.m22*b.y + A.m32*b.z,
		A.m13*b.x + A.m23*b.y + A.m33*b.z);
}

inline VECTOR3 operator* (const VECTOR3 &a, const double f)
{
	VECTOR3 c;
	c.x = a.x*f;
	c.y = a.y*f;
	c.z = a.z*f;
	return c;
}

inline VECTOR3 operator- (const VECTOR3 &a)
{
	VECTOR3 c;
	c.x = -a.x;
	c.y = -a.y;
	c.z = -a.z;
	return c;
}

inline VECTOR3 operator- (const VECTOR3 &a, const VECTOR3 &b)
{
	VECTOR3 c;
	c.x = a.x - b.x;
	c.y = a.y - b.y;
	c.z = a.z - b.z;
	return c;
}

inline double dotp(const VECTOR3 &a, const VECTOR3 &b)
{
	return a.x*b.x + a.y*b.y + a.z*b.z;
}

inline VECTOR3 crossp(const VECTOR3 &a, const VECTOR3 &b)
{
	return _V(a.y*b.z - b.y*a.z, a.z*b.x - b.z*a.x, a.x*b.y - b.x*a.y);
}

inline double length(const VECTOR3 &a)
{
	return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}

inline VECTOR3 operator/ (const VECTOR3 &a, const double f)
{
	VECTOR3 c;
	c.x = a.x / f;
	c.y = a.y / f;
	c.z = a.z / f;
	return c;
}

inline VECTOR3 unit(const VECTOR3 &a)
{
	return a / length(a);
}

inline VECTOR3 operator+ (const VECTOR3 &a, const VECTOR3 &b)
{
	VECTOR3 c;
	c.x = a.x + b.x;
	c.y = a.y + b.y;
	c.z = a.z + b.z;
	return c;
}

//Not martins code

inline MATRIX3 tmat(MATRIX3 a)
{
	MATRIX3 b;

	b = _M(a.m11, a.m21, a.m31, a.m12, a.m22, a.m32, a.m13, a.m23, a.m33);
	return b;
}

inline void agcSwap(VECTOR3 *a)
{
	double b = a->z; a->z = a->y; a->y = b;
}

inline double acos2(double a)
{
	if (a > 1.0)
	{
		a = 1.0;
	}
	else if (a < -1.0)
	{
		a = -1.0;
	}
	return acos(a);
}

inline double angle(VECTOR3 a, VECTOR3 b)
{
	return acos2(dotp(unit(a), unit(b)));
}

inline MATRIX3 _MRx(double a)
{
	double ca = cos(a), sa = sin(a);
	return _M(1.0, 0, 0, 0, ca, sa, 0, -sa, ca);
}

inline MATRIX3 _MRy(double a)
{
	double ca = cos(a), sa = sin(a);
	return _M(ca, 0, -sa, 0, 1.0, 0, sa, 0, ca);
}

inline MATRIX3 _MRz(double a)
{
	double ca = cos(a), sa = sin(a);
	return _M(ca, sa, 0, -sa, ca, 0, 0, 0, 1.0);
}