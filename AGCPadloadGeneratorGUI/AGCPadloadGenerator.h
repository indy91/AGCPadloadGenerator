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

struct IMUBiasCompensationData
{
	IMUBiasCompensationData();

	//PIPA X bias correction, cm/s
	double PBIASX;
	//PIPA X scale factor correct, ppm
	double PIPASCFX;
	//PIPA Y bias correction, cm/s
	double PBIASY;
	//PIPA Y scale factor correct, ppm
	double PIPASCFY;
	//PIPA Z bias correction, cm/s
	double PBIASZ;
	//PIPA z scale factor correct, ppm
	double PIPASCFZ;
	//X Gyro bias correction, meru
	double NBDX;
	//Y Gyro bias correction, meru
	double NBDY;
	//Z Gyro bias correction, meru
	double NBDZ;
	//X Gyro input axis drift, meru/g
	double ADIAX;
	//Y Gyro input axis drift, meru/g
	double ADIAY;
	//Z Gyro input axis drift, meru/g
	double ADIAZ;
	//X Gyro spin axis drift, meru/g
	double ADSRAX;
	//Y Gyro spin axis drift, meru/g
	double ADSRAY;
	//Z Gyro spin axis drift, meru/g
	double ADSRAZ;
};

struct BlockIICMCData
{
	//All CMC

	//P23 W-Matrix initial position error, ft
	double WMIDPOS;
	//P23 W-Matrix initial velocity error, ft/s
	double WMIDVEL;

	//All CMC except Colossus 237
	double EMDOT = 0.0;

	//Comanche 45 to 108
	double EK1VAL = 0.0;
	//Comanche 55 to 108
	double EK2VAL = 0.0;
	double EK3VAL = 0.0;
	double FANG = 0.0;
	//Artemis 72
	double EIMP1SEC = 0.0;
	double EFIMP01 = 0.0;
	double EFIMP16 = 0.0;
	//Comanche 108 and Artemis 72
	double TRUNSF = 0.0; // P/REV/CS
	double SHAFTSF = 0.0; //

	//Skylark only
	double C12ALPHA = 0.0; //Used in Docked DAP, SEC2/DEG
	double ECP = 0.0; //Used in docked DAP, unitless
	double ECYW = 0.0; //Used in docked DAP, unitless
	double ALPHAP = 0.0; //Used in docked DAP, unitless
	double ALPHAYW = 0.0; //Used in docked DAP, unitless
	double KMJDCKD = 0.0; //Used in docked DAP, DEG/SEC2
	double KMJ1DCKD = 0.0; //Used in docked DAP, DEG/SEC2
	double KMJ2DCKD = 0.0; //Used in docked DAP, DEG/SEC2
	double JMDCKD = 0.0; //Used in docked DAP, SEC2/DEG
	double JM1DCKD = 0.0; //Used in docked DAP, SEC2/DEG
	double JM2DCKD = 0.0; //Used in docked DAP, SEC2/DEG
	int CH6FAIL = 0; //R3 of N87. Docked DAP channel 6 jet inhibit mask
	double DKRATE = 0.0; //R1 of N89, Docked DAP maneuver rate, DEG/SEC
	double TTPI = 0.0; //Time of TPI, seconds

	IMUBiasCompensationData IMUBiasCompensation;
};

struct LGCData
{
	IMUBiasCompensationData IMUBiasCompensation;
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
	//P-20, ft/s
	double RVARMIN;

	//If true then use R2 gravity model (only supported by Open Orbiter)
	bool R2Model = false;

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

	struct PIOSDataSet
	{
		std::string name;
		int epoch = 0;
		double w_E = 0.0;
		double B_0 = 0.0;
		double Omega_I0 = 0.0;
		double F_0 = 0.0;
		double B_dot = 0.0;
		double Omega_I_dot = 0.0;
		double F_dot = 0.0;
		double cosI = 0.0;
		double sinI = 0.0;
		double t0 = 0.0;
		bool AZ0Hardcoded = false; //Later AGC version have AZ0 hardcoded
		double AZ0 = 0.0;
		bool ExtendedLimit = false; //false = 481 days limit on ephemerides and PIOS, true = 963 days limit
	};

	enum AGCVersions
	{
		AGCVersionError,
		Solarium055,
		Sundisk282,
		Colossus237,
		Colossus249,
		Comanche045,
		Comanche055,
		Comanche067,
		Comanche072,
		Comanche108,
		Artemis072,
		Sunburst120,
		Sundance306,
		Luminary069,
		Luminary069R2,
		Luminary099,
		Luminary116,
		Luminary131,
		Luminary131R1,
		Luminary178,
		Luminary210,
		Zerlina56,
		Skylark048
	};

public:
	AGCPadloadGenerator();
	~AGCPadloadGenerator();

	void RunBlockI();
	int RunCMC();
	int RunLGC();
	void GenerateRopeConstants(int Year);

	//Padloaded related variables
	double LaunchMJD;
	std::string Pad;
	std::string RopeName;
	std::string PIOSDataSetName;
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
	BlockIICMCData CMCDATA;
	LGCData LGCDATA;

protected:
	void SaveEMEM(int address, int value);
	void WriteEMEM(int address, int value, bool cmc);
	void AGCCorrectionVectors(PIOSDataSet dataset, double mjd_launchday, double dt_UNITW, double dt_504LM, bool IsCMC, bool IsSundance);
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

	AGCVersions GetCMCVersion(std::string name);
	AGCVersions GetLGCVersion(std::string name);

	int GetPIOSDataSet(std::string name, PIOSDataSet &data);

	std::vector<EMEM> arr;
	std::ofstream myfile, debugfile;
	char Buffer[256];

	double PrelaunchMJD;

	int iTemp, iTemp2, iTemp3;
	double dTemp;

	void SetPadData(Launchpad pad);

	void IMUCompensation(bool cmc, bool earlymodel);
	void R2Model(int address);

	//Same addresses for all CMCs
	void CMCDefaults(bool EarlyPIPABias, bool IsC108 = false);

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
	void LGCDefaults(bool EarlyPIPABias, bool mass = false);

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

	void Skylark048Padload();
	void SkylarkSolarEphemeris(double TC, double T0);
	void SkylarkCorrectionMatrix(double TC, double T0);

	//Block I
	void BlockIDefaults();
	void CoronaDefaults();
	void SolariumDefaults();

	//Conversions
	MATRIX3 J2000EclToBRCSMJD(double mjd) const;
	MATRIX3 J2000EclToBRCS(int epoch) const;
	double TJUDAT(int Y, int M, int D) const;
	MATRIX3 GetRotationMatrix(int plan, double t) const;
	double MJDOfNBYEpoch(int epoch) const;
	double MJD2JD(double MJD) const;
	double JD2MJD(double JD) const;

	Earth earth;
};