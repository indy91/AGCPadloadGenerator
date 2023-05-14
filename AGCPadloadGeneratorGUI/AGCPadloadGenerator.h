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

struct BlockIData
{
	//Interval between the launch vector passed through the inertial z-x plane and the time the AGC clock was zeroed at midnight of launch date
	double DTEPOCH;
	//Atlantic abourt site
	double T_ATL, lat_ATL, lng_ATL;
	//Pacific abort site
	double T_PAC, lat_PAC, lng_PAC;

	//SPS-1
	double e_SPS1, a_SPS1;
	//SPS-2
	double e_SPS2, a_SPS2;

	//Time from lift-off at which roll monitor begins
	double TROLL;
	//Time from lift-off at which pitch monitor begins
	double TPITCH;
	//Time pitch monitor is on
	double TENDPITCH;
	//Pitch polynomial
	double POLYCOFF[7];
};

class AGCPadloadGenerator
{
	enum Launchpad
	{
		LC39A,
		LC39B,
		LC34,
		LC37B
	};

public:
	AGCPadloadGenerator();
	~AGCPadloadGenerator();

	void RunBlockI();
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
	double PadLat, PadLong, PadAlt, PadAzimuth;
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
	//P20 W-Matrix maximum position deviation that gets processed automatically, m
	double RMAX;
	//P20 W-Matrix maximum velocity deviation that gets processed automatically, m
	double VMAX;

	BlockIData BLOCKI;

protected:
	void SaveEMEM(int address, int value);
	void WriteEMEM(int address, int value, bool cmc);
	void AGCCorrectionVectors(std::string rope, double mjd_launchday, double dt_UNITW, double dt_504LM, bool IsCMC);
	void AGCEphemeris(double T0, int Epoch, double TEphem0, double Span);

	int clbkEphemeris(int body, double mjd, int req, double *ret);
	int agcCelBody_RH(int Cel, double mjd, int Flags, VECTOR3 *Pos = NULL, VECTOR3 *Vel = NULL);
	int agcCelBody_LH(int Cel, double mjd, int Flags, VECTOR3 *Pos = NULL, VECTOR3 *Vel = NULL);

	MATRIX3 CalculateEarthTransformationMatrix(double t_M, double A_Z0, double w_E);
	MATRIX3 CalculateMoonTransformationMatrix(double t_M, double B_0, double B_dot, double Omega_I0, double Omega_I_dot, double F_0, double F_dot, double cosI, double sinI);
	VECTOR3 r_from_latlong(double lat, double lng);
	VECTOR3 r_from_latlong(double lat, double lng, double alt, int P, int F);

	MATRIX3 SolariumEarthFixedToSM(double lat, double lng, double azi);
	double Solarium055DTEPOCHCalculation(double A_Z0, double MJD_0, double MJD_L, double lng);
	double HANGLE(int E, int Y, int D);

	std::vector<EMEM> arr;
	std::ofstream myfile, debugfile;
	char Buffer[256];

	double PrelaunchMJD;

	int iTemp, iTemp2, iTemp3;
	double dTemp;

	void SetPadData(Launchpad pad);

	//Same addresses for all CMCs
	void CMCDefaults();

	//Colossus
	void Colossus237_249_Defaults(bool Is249);
	//Comanche
	void Comanche55Defaults();
	void Comanche67Defaults();
	//Artemis
	void ArtemisDefaults();

	//Same addresses for all LGCs
	void LGCDefaults();
	void Luminary131Defaults();

	//
	void BlockIDefaults();
	void CoronaDefaults();
	void SolariumDefaults();

	Earth earth;
};