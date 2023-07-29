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

struct BlockCMCIIData
{

};

struct BlockLGCIIData
{

};

struct BlockIIData
{
	double CSMMass; //lbs
	double LMMass;  //lbs
	double TotalMass; //lbs

	//P20 W-Matrix initial position error, ft
	double WRENDPOS;
	//P20 W-Matrix initial velocity error, ft/s
	double WRENDVEL;
	//P20 W-Matrix maximum position deviation that gets processed automatically, ft
	double RMAX;
	//P20 W-Matrix maximum velocity deviation that gets processed automatically, ft/s
	double VMAX;

	//CMC Only

	//Entry target data for boost aborts, latitude, degrees
	double LAT_SPL;
	//Entry target data for boost aborts, longitude, degrees
	double LNG_SPL;
	double PACTOFF;
	double YACTOFF;
	double LADPAD;
	double LODPAD;
	double ALFAPAD;
	double P37RANGE;
	//Pitch polynomial
	double POLYNUM[7];
	//Start of pitch polynomial in seconds since liftoff
	double RPSTART;
	//End of pitch polynomial in seconds since liftoff
	double POLYSTOP;

	//LGC Only
	double WSHAFT;
	double WTRUN;
	double SHAFTVAR;
	double TRUNVAR;
	double WSURFPOS;
	double WSURFVEL;
	double HIASCENT;
	double AGSK;
	double ROLLTIME;
	double PITCHTIME;
	double TLAND;

	//Descent constants
	VECTOR3 RBRFG;
	VECTOR3 VBRFG;
	VECTOR3 ABRFG;
	double VBRFG_star;
	double ABRFG_star;
	double JBRFG_star;
	double GAINBRAK;
	double TCGFBRAK;
	double TCGIBRAK;
	VECTOR3 RAPFG;
	VECTOR3 VAPFG;
	VECTOR3 AAPFG;
	double VAPFG_star;
	double AAPFG_star;
	double JAPFG_star;
	double GAINAPPR;
	double TCGFAPPR;
	double TCGIAPPR;
	double VIGN;
	double RIGNX;
	double RIGNZ;
	double KIGNXB4;
	double KIGNYB8;
	double KIGNVB4;
	double HIGHCRIT;
	double LOWCRIT;
	//Expected angular acceleration at lunar liftoff
	double IGNAOSQ;
	double IGNAOSR;
	//Terrain model
	double ABSC[5], SLOPE[5];
	//Descent abort
	double J1PARM, K1PARM, J2PARM, K2PARM, THETCRIT, RAMIN;
	//LR data reasonability test parameter, feet
	double DELQFIX;
	//An augment added to TTT during initialization for the first guidance cycle of the approach phase, seconds
	double DELTTFAP;
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

	BlockIData BLOCKI;
	BlockIIData BLOCKII;

protected:
	void SaveEMEM(int address, int value);
	void WriteEMEM(int address, int value, bool cmc);
	void AGCCorrectionVectors(std::string rope, double mjd_launchday, double dt_UNITW, double dt_504LM, bool IsCMC, bool IsSundance);
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
	void CMCDefaults(bool IsC108 = false);

	//Colossus
	void Colossus237_249_Defaults(bool Is249);
	//Comanche
	void Comanche45Padload(bool IsR2);
	void Comanche55Padload();
	void Comanche67Padload();
	void Comanche72Padload();
	void Comanche108Padload();
	//Artemis
	void Artemis72Padload();
	void SavePOLYNUM(int address);

	//Same addresses for all LGCs
	void LGCDefaults(bool mass = false);

	void Sundance306Defaults();
	void Luminary069Padload(bool IsR2);
	void Luminary099Padload();
	void Luminary116Padload();
	void Luminary131Padload();
	void Luminary178Padload();
	void Zerlina56Padload();
	void Luminary210Padload();

	void Luminary099_116_Defaults();
	void DescentConstants11_13();
	void DescentConstants14_17();

	//Block I
	void BlockIDefaults();
	void CoronaDefaults();
	void SolariumDefaults();

	Earth earth;
};