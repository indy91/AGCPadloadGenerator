// Copyright (c) Martin Schweiger
// Licensed under the MIT License

#include "Vsop87.h"
#include <stdio.h>
#include <ctype.h>
#include <fstream>
#include <math.h>

#define DLLCLBK extern "C" __declspec(dllexport)

using namespace std;

// ===========================================================
// class VSOPOBJ
// Base class for planets controlled by VSOP87 solutions
// ===========================================================

VSOPOBJ::VSOPOBJ ()
{
	a0 = 1.0;               // should be overwritten by derived class
	prec = 1e-6;            // default precision
	termidx = 0;
	termlen = 0;
	term = 0;
	//sp[0].t = sp[1].t = -1e20; // invalidate
	SetSeries ('B');        // default series: spherical, J2000
}

VSOPOBJ::~VSOPOBJ ()
{
	if (termidx) delete []termidx;
	if (termlen) delete []termlen;
	if (term) delete []term;
}

void VSOPOBJ::clbkInit (double ErrorLimit)
{
	prec = ErrorLimit;
}

void VSOPOBJ::SetSeries (char series)
{
	sid = toupper (series);
	fmtflag = 0;
	if (sid == 'B' || sid == 'D') fmtflag |= EPHEM_POLAR;
	if (sid == 'E')               fmtflag |= EPHEM_PARENTBARY;
}

bool VSOPOBJ::ReadData (const char *name)
{
	int nterm, cooidx, alpha, i, iused, nused = 0, ntot = 0;
	double a, b, c, tfac, err;
	TERM3 **ppterm, *pterm;

	char cbuf[256];
	sprintf_s(cbuf, "Vsop87%c.dat", sid);
	ifstream ifs (cbuf);
	if (!ifs) {
		//oapiWriteLogError("VSOP87 %s: Data file not found: %s", name, cbuf);
		return false;
	}

	ifs >> nalpha;

	ppterm  = new TERM3*[(nalpha+1)*3];
	termidx = new IDX3[nalpha+1];
	termlen = new IDX3[nalpha+2];

	for (cooidx = 0; cooidx < 3; cooidx++) {
		tfac = 1.0;
		for (alpha = 0; alpha <= nalpha; alpha++) {
			ifs >> nterm;
			pterm = ppterm[cooidx*(nalpha+1)+alpha] = new TERM3[nterm];
			for (i = 0, iused = nterm; i < nterm; i++) {
				ifs >> a >> b >> c;
				if (iused == nterm) {
					pterm[i][0] = a, pterm[i][1] = b, pterm[i][2] = c;
					if (cooidx == 2) a /= a0; // radius in terms of mean SMa
					err = 2.0*sqrt (i+1.0)*a*tfac;
					if (err < prec) iused = i;
				}
			}
			termlen[alpha][cooidx] = iused;
			termidx[alpha][cooidx] = nused;
			nused += iused;
			ntot  += nterm;
			tfac *= 5.0; // don't ask
		}
		termlen[alpha][cooidx] = 0;
	}
	// now copy everything into a single array
	term = new TERM3[nused];
	for (cooidx = 0; cooidx < 3; cooidx++) {
		for (alpha = 0; alpha <= nalpha; alpha++) {
			memcpy (term+termidx[alpha][cooidx], ppterm[cooidx*(nalpha+1)+alpha], termlen[alpha][cooidx]*sizeof(TERM3));
			delete []ppterm[cooidx*(nalpha+1)+alpha];
		}
	}
	delete []ppterm;

	//oapiWriteLogV("VSOP87(%c) %s: Precision %0.1le, Terms %d/%d", sid, name, prec, nused, ntot);
	return true;
}

// ===========================================================
// Name: VsopEphem()
// Desc: Return heliocentric ecliptic spherical positions and
//       velocity for planet 'obj' (VSOP_MERCURY to VSOP_NEPTUNE
//       at time 'mjd' for ecliptic and equinox J2000.
//       Values returned in 'ret' are
//       ret[0] = longitude l [rad]
//       ret[1] = latitude b  [rad]
//       ret[2] = radius r    [AU]
//       ret[3] = velocity in longitude [rad/s]
//       ret[4] = velocity in latitude  [rad/s]
//       ret[5] = radial velocity [AU/s]
// ===========================================================
void VSOPOBJ::VsopEphem (double mjd, double *ret)
{
	static const double mjd2000 = 51544.5;  // MJD date of epoch J2000
	static const double a1000   = 365250.0; // days per millenium
	static const double rsec    = 1.0/(a1000*86400.0); // 1/seconds per millenium

	static const double c0   = 299792458;     // speed of light [m/s]
	static const double tauA = 499.004783806; // light time for 1 AU [s]
	static const double AU   = c0*tauA;       // 1 AU in meters
	static const double pscl = AU;            // convert AU -> m
	static const double vscl = AU*rsec;       // convert AU/millenium -> m/s

	TERM3 *pterm;

	double a, b, c, arg, tm, termdot;
	int i, cooidx, alpha;

	// zero result array
	for (i = 0; i < 6; i++) ret[i] = 0.0;

	// set time and powers
	double t[VSOP_MAXALPHA+1];
	t[0] = 1.0;
	t[1] = (mjd-mjd2000)/a1000;
	for (i = 2; i <= VSOP_MAXALPHA; ++i) t[i] = t[i-1] * t[1];

	// term summation
	for (cooidx = 0; cooidx < 3; ++cooidx) { // loop over spatial dimensions

		for (alpha = 0; termlen[alpha][cooidx]; ++alpha) { // loop over powers of time

			pterm = term+termidx[alpha][cooidx];
			tm = termdot = 0.0;
		    for (i = 0; i < termlen[alpha][cooidx]; ++i) {
				a        = pterm[i][0];
				b        = pterm[i][1];
				c        = pterm[i][2];
				arg      = b + c * t[1];
				tm      += a * cos(arg);
				termdot -= c * a * sin(arg);
			}
			ret[cooidx] += t[alpha] * tm;
			ret[cooidx+3] += t[alpha] * termdot +
				(alpha > 0 ? alpha * t[alpha - 1] * tm : 0.0);

		} // end loop alpha
	} // end loop cooidx

	if (fmtflag & EPHEM_POLAR) {
		// convert millenium rate to second rate
		for (i = 3; i < 6; i++) ret[i] *= rsec;
		// should also convert radius from AU to m
	} else {
		double tmp;
		for (i = 0; i < 3; i++) ret[i] *= pscl;
		for (     ; i < 6; i++) ret[i] *= vscl;
		// swap y and z to map to orbiter system
		tmp = ret[1]; ret[1] = ret[2]; ret[2] = tmp;
		tmp = ret[4]; ret[4] = ret[5]; ret[5] = tmp;
	}
}

