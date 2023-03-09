// Copyright (c) Martin Schweiger
// Licensed under the MIT License

#define ORBITER_MODULE

#include "Earth.h"

// ======================================================================
// class Earth: implementation
// ======================================================================

Earth::Earth (): VSOPOBJ ()
{
	a0         = 1.0;   // semi-major axis [AU]
}

void Earth::clbkInit ()
{
	VSOPOBJ::clbkInit(1.e-8);
	ReadData ("Earth");
}

int Earth::clbkEphemeris (double mjd, int req, double *ret)
{
	VsopEphem (mjd, ret);
	return fmtflag | EPHEM_TRUEPOS | EPHEM_TRUEVEL;
}