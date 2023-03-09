#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include "VectorMath.h"
#include "Earth.h"

struct EMEM
{
	int address;
	int value;

	//For sorting
	bool operator<(const EMEM& rhs) const
	{
		return address < rhs.address;
	}
};

class AGCPadloadGenerator
{
public:
	AGCPadloadGenerator();
	~AGCPadloadGenerator();

	void RunCMC();
	void RunLGC();

	//Padloaded related variables
	double LaunchMJD;
	std::string Pad, RopeName;
	double T_504LM, T_UNITW, EphemerisSpan;
	double LSLat, LSLng, LSAlt;
	double TLAND;
	double RTED1, EMSALT;
	double LaunchAzimuth;
	//Variance of the primary body radius vector, m^2
	double RPVAR;
	//Initial condition for diagonal elements of W8 for known landmark, m
	double S22WSUBL;
	//Initial position error for landmark tracking, m
	double WORBPOS;
	//Initial velocity error for landmark tracking, m/s
	double WORBVEL;
	double PadLat, PadLong, PadAlt;
	//Optical alignment target
	double TAZEL[4];
	double TEPHEM;
	VECTOR3 RLS;
	//Limit for marks, 3 bits of CDU angle changes in this time, s
	double CDUCHKWD;
	//Threshold value of DV which must be sensed in two second period or SPS thrust failure is indicated, m/s
	double DVTHRESH;
	//Horizon altitude for P23, m
	double HORIZALT;
	//Alternate line-of-sight variance, P20, rad^2
	double ALTVAR;
	//P20 W-Matrix initial position error, m
	double WRENDPOS;
	//P20 W-Matrix initial velocity error, m
	double WRENDVEL;
protected:
	void SaveEMEM(int address, int value);
	void WriteEMEM(int address, int value, bool cmc);
	void AGCCorrectionVectors(std::string rope, double mjd_launchday, double dt_UNITW, double dt_504LM, bool IsCMC);
	void AGCEphemeris(double T0, int Epoch, double TEphem0, double Span);

	int clbkEphemeris(int body, double mjd, int req, double *ret);
	int agcCelBody_RH(int Cel, double mjd, int Flags, VECTOR3 *Pos = NULL, VECTOR3 *Vel = NULL);
	int agcCelBody_LH(int Cel, double mjd, int Flags, VECTOR3 *Pos = NULL, VECTOR3 *Vel = NULL);

	VECTOR3 CalculateRLS(double lat, double lng, double alt, double rad);
	MATRIX3 CalculateEarthTransformationMatrix(double t_M, double A_Z0, double w_E);
	MATRIX3 CalculateMoonTransformationMatrix(double t_M, double B_0, double B_dot, double Omega_I0, double Omega_I_dot, double F_0, double F_dot, double cosI, double sinI);

	std::vector<EMEM> arr;
	std::ofstream myfile, debugfile;
	char Buffer[256];

	double PrelaunchMJD;

	int iTemp, iTemp2, iTemp3;
	double dTemp;

	//Same addresses for all CMCs
	void CMCDefaults();

	//Colossus
	void Colossus237Defaults();
	//Comanche
	void Comanche55Defaults();
	void Comanche67Defaults();
	//Artemis
	void ArtemisDefaults();

	//Same addresses for all LGCs
	void LGCDefaults();
	void Luminary131Defaults();


	Earth earth;
};