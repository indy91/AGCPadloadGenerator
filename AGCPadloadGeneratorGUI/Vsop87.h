// Copyright (c) Martin Schweiger
// Licensed under the MIT License

#pragma once

#define VSOP_MAXALPHA 5		// max power of time

#define EPHEM_TRUEPOS     0x01	///< true body position
#define EPHEM_TRUEVEL     0x02	///< true body velocity
#define EPHEM_BARYPOS     0x04	///< barycentric position
#define EPHEM_BARYVEL     0x08	///< barycentric velocity
#define EPHEM_BARYISTRUE  0x10	///< body has no child objects
#define EPHEM_PARENTBARY  0x20	///< ephemerides are computed in terms of the barycentre of the parent body's system
#define EPHEM_POLAR       0x40	///< data is returned in polar format

typedef int IDX3[3];
typedef double TERM3[3];

// ===========================================================
// class VSOPOBJ
// Base class for planets controlled by VSOP87 solutions
// ===========================================================

class VSOPOBJ {
public:
	VSOPOBJ ();
	virtual ~VSOPOBJ ();
	void clbkInit (double ErrorLimit);

protected:
	void SetSeries (char series);
	// Set VSOP series ('A' to 'E')

	bool ReadData (const char *name);
	// Read perturbation terms up to required accuracy from data file

	void VsopEphem (double mjd, double *ret);
	// Calculate ephemerides


	double a0;       // semi-major axis [AU]
	double prec;     // tolerance limit (1e-3 .. 1e-8)
	int fmtflag;     // data format flag
	int nalpha;      // order of time polynomials
	IDX3 *termidx;   // term index list
	IDX3 *termlen;   // term list lengths
	TERM3 *term;     // term list
	//Sample sp[2];

private:
	//void Interpolate (double t, double *data, const Sample *s0, const Sample *s1);

	char sid;
	int datatp;  // return data type: true pos or barycentric
};