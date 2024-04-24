#include "AGCPadloadGenerator.h"
#include "VectorMath.h"

VECTOR3 StarCatalogue[] =
{
	_V(0.87325707, 0.222717753, 0.433380771),
	_V(0.933983515, 0.0421048982, -0.354826677),
	_V(0.474230235, 0.456854026, 0.75258892),
	_V(0.492731712, -0.129171941, -0.860540568),
	_V(0.010128911, 0.40497893, 0.914269912),
	_V(0.543321089, 0.233682862, -0.806346398),
	_V(0.698256178, 0.681882483, -0.217886645),
	_V(0.404981298, 0.764255442, 0.501900156),
	_V(0.343908924, 0.93415959, -0.0952497353),
	_V(0.195051858, 0.833544498, -0.516873623),
	_V(0.130502244, 0.912114029, 0.388609268),
	_V(-0.0632195158, 0.236597063, -0.969548928),
	_V(-0.187545554, 0.747401122, -0.637352514),
	_V(-0.418201486, 0.865470289, -0.27580554),
	_V(-0.362955763, 0.232568134, -0.902316561),
	_V(-0.471158209, 0.73100027, 0.493607686),
	_V(-0.777933252, 0.49981843, -0.380790484),
	_V(-0.864519888, 0.502532948, 0.00812396612),
	_V(-0.96672924, 0.142390103, 0.212507964),
	_V(-0.951226902, -0.180199054, -0.250391056),
	_V(-0.449404749, -0.402805961, -0.797359849),
	_V(-0.914079394, -0.403947962, -0.0358455889),
	_V(-0.581451821, 0.0311390298, 0.812984712),
	_V(-0.685292494, -0.623812256, -0.375809082),
	_V(-0.783614901, -0.352762709, 0.511376728),
	_V(-0.529185544, -0.481419898, 0.698711344),
	_V(-0.344817825, -0.935282024, -0.0796756118),
	_V(-0.109619054, -0.684031025, -0.721169342),
	_V(-0.106559622, -0.803652501, 0.585480746),
	_V(0.125067322, -0.456806946, 0.88073014),
	_V(0.214089696, -0.974959511, -0.0601627197),
	_V(0.459149324, -0.741396708, 0.489400469),
	_V(0.558073669, -0.825925685, 0.0800033956),
	_V(0.325589471, -0.737606353, -0.591547432),
	_V(0.455646022, -0.209297007, 0.865206025),
	_V(0.817196194, -0.436629484, 0.376224766),
	_V(0.837327731, -0.410968743, -0.360537048)
};

void FormatStarAGCUnitVector(char *Buffer, double val)
{
	//Gets rid of the leading zero

	char Buffer2[128];

	sprintf_s(Buffer2, "%+.10lf", val);

	int i, j;

	j = 0;

	for (i = 0; i < 14; i++)
	{
		if (i != 1)
		{
			Buffer[j] = Buffer2[i];
			j++;
		}
	}
}

void FormatAGCUnitVector(char *Buffer, double val)
{
	//Gets rid of the leading zero

	char Buffer2[128];

	sprintf_s(Buffer2, "%.9lf", val);

	int i, j, num;

	j = 0;

	if (val >= 0.0)
	{
		num = 0;
	}
	else
	{
		num = 1;
	}

	for (i = 0; i < 14; i++)
	{
		if (i != num)
		{
			Buffer[j] = Buffer2[i];
			j++;
		}
	}
}

void FormatAGCConstants(double val, double &val2, int &num)
{
	num = 0;
	val2 = val;

	while (abs(val2) < 1.0)
	{
		num--;
		val2 = val * pow(10.0, -num);
	}
}

double AngleBetweenMatrices(MATRIX3 Rot1, MATRIX3 Rot2)
{
	MATRIX3 Rot3;
	double Tr;

	Rot3 = mul(tmat(Rot1), Rot2);
	Tr = Rot3.m11 + Rot3.m22 + Rot3.m33; //Trace

	return acos((Tr - 1.0) / 2.0);
}

void CalculateSolarEphemerisConstants(double TC, double T0, double &K1, double &K3, double &LOSO, double &LOSR, double &CARG, double &CMOD, double &OMEGAC)
{
	//From SGA Memo #8-71 and Addendum

	//TC: julian ephemeris date at midpoint of the time interval over which the solar position approximation is desired
	//T0: julian ephemeris date midnight July 1 before launch

	double D, fact, n_EM, eps, e_EM, W_0EM, M_0EM, fact2;

	//Time variables
	D = (TC - 2415020.0) / 1e4;
	fact = 1.0 + 2.852e-8; //Converts ephemeris time rate to universal time rate
	n_EM = fact * 0.985600267; //VAL67 + 8 (1/365)

	fact2 = 0.32328 / 36525.0; //Precession correction in degrees per ephemeris day

	//Calculate mean obliquity
	eps = 23.452294 - 3.5626e-3*D - 1.23e-7*D*D + 1.03e-8*D*D*D;
	eps *= RAD;
	//Calculate mean eccentricity
	e_EM = 0.01675104 - 1.1444e-5*D - 9.4e-9*D*D;
	//Calculate longitude of perihelion
	W_0EM = 101.220833 + 0.470684*D - fact2*(TC - T0) + 3.39e-5*D*D + 7.0e-8*D*D*D;
	W_0EM = fmod(W_0EM, 360.0);
	W_0EM *= RAD;

	//Calculate mean anomaly
	M_0EM = 358.475845 + n_EM * (T0 - 2415020.0) - 1.12e-5*D*D - 7.0e-8*D*D*D;
	M_0EM = fmod(M_0EM, 360.0);
	M_0EM *= RAD;

	K1 = cos(eps); //KONMAT+4
	K3 = sin(eps); //KONMAT+6

	LOSO = W_0EM + M_0EM - PI; //LOSO
	LOSR = (n_EM + fact2)*RAD;

	CMOD = (2.0 * e_EM - 0.25*pow(e_EM, 3)); //C
	OMEGAC = n_EM*RAD;
	CARG = M_0EM - PI; // PHASEC
}

double MoonMeanLongitudeAscendingNode(double T0)
{
	//T0 = julian ephemeris date

	double T, d, D;

	T = (T0 - 2415020.0) / 36525.0;
	D = 3.6525*T;
	d = 10000.0*D;

	double Omega;
	Omega = 259.183275 - 0.0529539222*d + 0.0001557*pow(D, 2) + 0.00000005*(D, 3);

	Omega = fmod(Omega, 360.0);
	if (Omega < 0.0)
	{
		Omega += 360.0;
	}

	Omega *= RAD;

	return Omega;
}

void CalculateLunarEphemerisConstants(double TC, double T0, double K1, double K3, double &K2, double &K4, double &LOMO, double &LOMR, double &A, double &B, double &OMEGAA, double &PHASEA, double &OMEGAB, double &PHASEB, double &LONO, double &LONR)
{
	//TC: julian ephemeris date at midpoint of the time interval over which the solar position approximation is desired
	//T0: julian ephemeris date midnight July 1 before launch

	double IM, LOMO_1972, PHASEA_1972, PHASEB_1972, JD_1972, DT;

	IM = 5.145396374132591*RAD; //Constant in 1961 explanatory supplement to the Astronomical Ephemeris

	K2 = -K3 * sin(IM);
	K4 = K1 * sin(IM);

	OMEGAA = 0.036291713*PI2;		// rad/day, sidereal period of the moon
	OMEGAB = 0.03125*PI2;			// rad/day, second term in Brown's series for the position of the Moon

	//TBD: Calculate lunar ephemeris constants (below the 1972 values)
	LOMO_1972 = 0.534104635*PI2;	// rad
	LOMR = 0.036600997*PI2;			// rad/day
	A = 0.017519236*PI2;			// rad
	PHASEA_1972 = 0.02452343*PI2;	// rad
	B = 0.003473642*PI2;			// rad
	PHASEB_1972 = 0.540564756*PI2;	// rad
	JD_1972 = 2441133.5;

	//Adjust for new epoch
	DT = T0 - JD_1972;

	LOMO = LOMO_1972 + LOMR * DT;
	LOMO = fmod(LOMO, PI2);

	LONO = MoonMeanLongitudeAscendingNode(T0);
	LONR = (MoonMeanLongitudeAscendingNode(T0 + 370.0) - LONO) / 370.0;

	PHASEA = PHASEA_1972 + OMEGAA * DT;
	PHASEA = fmod(PHASEA, PI2);

	PHASEB = PHASEB_1972 + OMEGAB * DT;
	PHASEB = fmod(PHASEB, PI2);
}

VECTOR3 EvaluateSolarEphemeris(double t, double LOSO, double LOSR, double C, double OMEGAC, double PHASEC, double K1, double K3)
{
	VECTOR3 U_ES;
	double LOS;

	LOS = LOSO + LOSR * t - C * sin(OMEGAC*t + PHASEC);
	U_ES = _V(cos(LOS), K1*sin(LOS), K3*sin(LOS));

	return unit(U_ES);
}

VECTOR3 EvaluateLunarEphemeris(double t, double LOMO, double LOMR, double A, double OMEGAA, double PHASEA, double B, double OMEGAB, double PHASEB, double LONO, double LONR, double K1, double K2, double K3, double K4)
{
	VECTOR3 U_EM;
	double LOM, LON;

	LOM = LOMO + LOMR * t - (A*sin(OMEGAA*t + PHASEA) + B * sin(OMEGAB*t + PHASEB));
	LON = LONO + LONR * t;

	U_EM = _V(cos(LOM), K1*sin(LOM) + K2 * sin(LOM - LON), K3*sin(LOM) + K4 * sin(LOM - LON));

	return unit(U_EM);
}

void AGCPadloadGenerator::GenerateRopeConstants(int Year)
{
	MATRIX3 Rot;
	double MJD_C;
	char Buffer[256], Buffer2[256];

	//MJD at epoch
	MJD_C = MJDOfNBYEpoch(Year);
	//Rotation matrix from J2000 ecliptic to BRCS
	Rot = J2000EclToBRCSMJD(MJD_C);

	std::ofstream file;
	VECTOR3 vTemp;

	file.open("RopeConstants.txt");

	file << "AGC CONSTANTS FOR THE YEAR " << Year << std::endl << std::endl;

	double valtemp;
	int inttemp, i, j;

	//STAR TABLES
	file << "Star Tables for CMC and LGC:" << std::endl;

	for (i = 36; i >= 0; i--)
	{
		vTemp = mul(Rot, StarCatalogue[i]);

		for (j = 0; j < 3; j++)
		{
			FormatStarAGCUnitVector(Buffer, vTemp.data[j]);
			file << "2DEC " << Buffer << " B-1 #STAR " << i + 1 << " ";
			if (j == 0)
			{
				file << "X";
			}
			else if (j == 1)
			{
				file << "Y";
			}
			else
			{
				file << "Z";
			}
			file << std::endl;
		}
		file << std::endl;
	}

	//TBD: EARTH PIOS
	file << "Earth Planetary Inertial Orientation Subroutine Constants:" << std::endl;

	MATRIX3 J2000, J2000_2, R, R3, R3_2;
	double MJD_0, A_Z, A_Z0, w_E;

	//TBD: Make variable?
	w_E = 7.29211514667e-5;

	FormatAGCConstants(w_E / (100.0*PI2), valtemp, inttemp);
	sprintf_s(Buffer, 255, "WEARTH 2DEC* %.9lf E%d B23* # REVS/CSEC", valtemp, inttemp);
	file << Buffer << std::endl;

	//MJD at midnight July 1st before epoch
	MJD_0 = JD2MJD(TJUDAT(Year - 1, 7, 1));
	//J2000 = matrix converting from J2000 to BRCS
	J2000 = J2000EclToBRCS(Year);
	//J2000 = matrix converting from J2000 to TC epoch
	J2000_2 = J2000EclToBRCSMJD(MJD_C);
	//R = matrix converting from PCS to J2000
	R = MatrixRH_LH(GetRotationMatrix(BODY_EARTH, MJD_C));
	//R3 = PCS to BRCS
	R3 = mul(J2000, R);
	R3_2 = mul(J2000_2, R);
	A_Z = atan2(R3_2.m21, R3_2.m11);
	A_Z0 = A_Z - w_E * (MJD_C - MJD_0) * 24.0 * 3600.0;
	while (A_Z0 < 0.0)
	{
		A_Z0 += PI2;
	}

	FormatAGCConstants(A_Z0 / PI2, valtemp, inttemp);
	sprintf_s(Buffer, 255, "AZO 2DEC* %.9lf E%d B 0* # REVS", valtemp, inttemp);
	file << Buffer << std::endl << std::endl;

	//MOON PIOS
	file << "Moon Planetary Inertial Orientation Subroutine Constants:" << std::endl;

	PIOSDataSet dataset;

	dataset.w_E = w_E;
	dataset.AZ0 = A_Z0;

	//Assume Artemis and Luminary 1E
	dataset.AZ0Hardcoded = true;
	dataset.ExtendedLimit = true;

	//Use 1972 data set as reference
	double MJD_1972, B_0_1972, B_dot_1972, Omega_I0_1972, Omega_I_dot_1972, F_0_1972, F_dot_1972, I_1972;

	MJD_1972 = 41133.0;
	B_0_1972 = 4.09157363336e-1;
	B_dot_1972 = -7.19758599677e-14;
	Omega_I0_1972 = 5.521857147;
	Omega_I_dot_1972 = -1.070470131e-8;
	F_0_1972 = 4.11720655556;
	F_dot_1972 = 2.67240425480e-6;
	I_1972 = 5521.5 / 3600.0 * RAD;

	//Analytical model to generate Moon PIOS constants
	double I, DT, Precession_Rate;

	//Angle between the mean lunar equatorial plane and the plane of the ecliptic. Assume it stays constant.
	I = I_1972;

	dataset.cosI = cos(I);
	dataset.sinI = sin(I);

	//Time in seconds between July 1st before epoch and MJD when 1972 data set is valid
	DT = (MJD_0 - MJD_1972)*24.0*3600.0;

	//Calculate F_0 and B_0 as linear projections from 1972 data set
	dataset.B_0 = B_0_1972 + B_dot_1972 * DT;
	dataset.B_0 = fmod(dataset.B_0, PI2);

	dataset.F_0 = F_0_1972 + F_dot_1972 * DT;
	dataset.F_0 = fmod(dataset.F_0, PI2);

	//Do the same for Omega_I, except it needs to account for precession of the coordinate system
	Precession_Rate = PI2 / (9413040.4 * 24.0 * 3600.0);
	dataset.Omega_I0 = Omega_I0_1972 + Omega_I_dot_1972*DT + Precession_Rate * DT;
	dataset.Omega_I0 = fmod(dataset.Omega_I0, PI2);

	//Assume the rates stay constant
	dataset.B_dot = B_dot_1972;
	dataset.Omega_I_dot = Omega_I_dot_1972;
	dataset.F_dot = F_dot_1972;

	//Write constants
	FormatAGCConstants(dataset.cosI, valtemp, inttemp);
	sprintf_s(Buffer, 255, "COSI 2DEC* %.9lf E%d B-1* # COS (5521.5 SEC.)", valtemp, inttemp);
	file << Buffer << std::endl;
	FormatAGCConstants(dataset.sinI, valtemp, inttemp);
	sprintf_s(Buffer, 255, "SINI 2DEC* %.9lf E%d B-1* # SIN (5521.5 SEC.)", valtemp, inttemp);
	file << Buffer << std::endl;
	FormatAGCConstants(dataset.Omega_I_dot / (100.0*PI2), valtemp, inttemp);
	sprintf_s(Buffer, 255, "NODDOT 2DEC* %.9lf E%d B28* # REV/CSEC", valtemp, inttemp);
	file << Buffer << std::endl;
	FormatAGCConstants(dataset.F_dot / (100.0*PI2), valtemp, inttemp);
	sprintf_s(Buffer, 255, "FDOT 2DEC* %.9lf E%d B27* # REV/CSEC", valtemp, inttemp);
	file << Buffer << std::endl;
	FormatAGCConstants(dataset.B_dot / (100.0*PI2), valtemp, inttemp);
	sprintf_s(Buffer, 255, "BDOT 2DEC* %.9lf E%d B28* # REV/CSEC", valtemp, inttemp);
	file << Buffer << std::endl;
	FormatAGCConstants(dataset.Omega_I0 / PI2, valtemp, inttemp);
	sprintf_s(Buffer, 255, "NODIO 2DEC* %.9lf E%d B0* # REV", valtemp, inttemp);
	file << Buffer << std::endl;
	FormatAGCConstants(dataset.F_0 / PI2, valtemp, inttemp);
	sprintf_s(Buffer, 255, "FSUBO 2DEC* %.9lf E%d B0* # REV", valtemp, inttemp);
	file << Buffer << std::endl;
	FormatAGCConstants(dataset.B_0 / PI2, valtemp, inttemp);
	sprintf_s(Buffer, 255, "BSUBO 2DEC* %.9lf E%d B0* # REV", valtemp, inttemp);
	file << Buffer << std::endl;

	//Debug
	MATRIX3 RM, Rot1, Rot2;
	double MJD, t_M, ang;

	MJD = MJD_0;

	file << std::endl;
	file << "Maximum angular and distance errors at selected MJDs:" << std::endl;

	while (MJD < MJD_0 + 370.0)
	{
		//Check rotations by comparing rotation matrices from BRCS to Moon fixed

		//Time since July 1st, midnight
		t_M = (MJD - MJD_0) * 24.0 * 3600.0;
		//Left-handed Moon-fixed to ecliptic at MJD
		RM = GetRotationMatrix(BODY_MOON, MJD);
		RM = MatrixRH_LH(RM);

		//Orbiter
		Rot1 = tmat(mul(J2000, RM));

		//AGC
		Rot2 = CalculateMoonTransformationMatrix(t_M, dataset.B_0, dataset.B_dot, dataset.Omega_I0, dataset.Omega_I_dot, dataset.F_0, dataset.F_dot, dataset.cosI, dataset.sinI);

		//Compare
		ang = AngleBetweenMatrices(Rot1, Rot2);

		sprintf_s(Buffer, 256, "%.3f = %.8lf (degrees), %.3lf (meters)\n", MJD, ang*DEG, ang*1.73809e6);
		file << Buffer;

		MJD += 10.0;
	}

	//PIOS DATA SET FOR CONFIG FILE
	file << std::endl << "PIOS data set for the PIOSDataSets.txt file:" << std::endl;
	sprintf_s(Buffer, 255, "NBY%d %d %.10e %.10lf %.10lf %.10lf %.10e %.10e %.10e %.10lf %.10lf %.1lf %d %.10lf %d", Year, Year, dataset.w_E, dataset.B_0, dataset.Omega_I0, dataset.F_0, dataset.B_dot, dataset.Omega_I_dot, dataset.F_dot,
		dataset.cosI, dataset.sinI, MJD_0, dataset.AZ0Hardcoded, dataset.AZ0, dataset.ExtendedLimit);
	file << Buffer << std::endl << std::endl;

	//LGC EPHEMERIDES
	file << "LGC Solar and Lunar Ephemerides" << std::endl;

	double TC, T0;
	double K1, K3, LOSO, LOSR, CARG, CMOD, OMEGAC;
	double K2, K4, LOMO, LOMR, A, OMEGAA, PHASEA, B, OMEGAB, PHASEB, LONO, LONR;

	//Convert from MJD to JD
	TC = MJD_C + 2400000.5;
	T0 = MJD_0 + 2400000.5;

	//Calculate constants
	CalculateSolarEphemerisConstants(TC, T0, K1, K3, LOSO, LOSR, CARG, CMOD, OMEGAC);
	CalculateLunarEphemerisConstants(TC, T0, K1, K3, K2, K4, LOMO, LOMR, A, B, OMEGAA, PHASEA, OMEGAB, PHASEB, LONO, LONR);

	sprintf_s(Buffer, 255, "KONMAT 2DEC  1.0         B-1  #        *************");
	file << Buffer << std::endl;
	sprintf_s(Buffer, 255, "       2DEC  0           B-28 #                    *");
	file << Buffer << std::endl;
	sprintf_s(Buffer, 255, "       2DEC  0           B-28 #                    *");
	file << Buffer << std::endl;
	sprintf_s(Buffer, 255, "       2DEC  0           B-28 #                    *");
	file << Buffer << std::endl;
	FormatAGCUnitVector(Buffer2, K1);
	sprintf_s(Buffer, 255, "       2DEC* %s B-1* #  K1 = COS(OBL)", Buffer2);
	file << Buffer << std::endl;
	FormatAGCUnitVector(Buffer2, K2);
	sprintf_s(Buffer, 255, "       2DEC* %s B-1* #  K2 = SIN (OBL) SIN (IM) (-1)", Buffer2);
	file << Buffer << std::endl;
	sprintf_s(Buffer, 255, "       2DEC  0           B-28 #                    *");
	file << Buffer << std::endl;
	FormatAGCUnitVector(Buffer2, K3);
	sprintf_s(Buffer, 255, "       2DEC* %s B-1* #  K3 = SIN(OBL)", Buffer2);
	file << Buffer << std::endl;
	FormatAGCUnitVector(Buffer2, K4);
	sprintf_s(Buffer, 255, "       2DEC* %s B-1* #  K4 = COS (OBL) SIN (IM)", Buffer2);
	file << Buffer << std::endl;
	file << std::endl;

	FormatAGCUnitVector(Buffer2, LOMR / PI2);
	sprintf_s(Buffer, 255, "RATESP 2DEC* %s B+4* # LOMR", Buffer2);
	file << Buffer << std::endl;

	FormatAGCUnitVector(Buffer2, LOSR / PI2);
	sprintf_s(Buffer, 255, "       2DEC* %s B+4* # LOSR", Buffer2);
	file << Buffer << std::endl;

	FormatAGCUnitVector(Buffer2, LONR / PI2);
	sprintf_s(Buffer, 255, "       2DEC* %s B+4* # LONR", Buffer2);
	file << Buffer << std::endl;

	FormatAGCUnitVector(Buffer2, LOMO / PI2);
	sprintf_s(Buffer, 255, "       2DEC  %s      # LOMO", Buffer2);
	file << Buffer << std::endl;

	FormatAGCUnitVector(Buffer2, LOSO / PI2);
	sprintf_s(Buffer, 255, "       2DEC  %s      # LOSO", Buffer2);
	file << Buffer << std::endl;

	FormatAGCUnitVector(Buffer2, LONO / PI2);
	sprintf_s(Buffer, 255, "       2DEC  %s      # LONO", Buffer2);
	file << Buffer << std::endl << std::endl;

	FormatAGCUnitVector(Buffer2, A / PI2);
	sprintf_s(Buffer, 255, "VAL67  2DEC* %s B+1* # AMOD", Buffer2);
	file << Buffer << std::endl;

	FormatAGCUnitVector(Buffer2, PHASEA / PI2);
	sprintf_s(Buffer, 255, "       2DEC  %s      # AARG", Buffer2);
	file << Buffer << std::endl;

	FormatAGCUnitVector(Buffer2, OMEGAA / PI2);
	sprintf_s(Buffer, 255, "       2DEC* %s B+1* # 1/27", Buffer2);
	file << Buffer << std::endl;

	FormatAGCUnitVector(Buffer2, B / PI2);
	sprintf_s(Buffer, 255, "       2DEC* %s B+1* # BMOD", Buffer2);
	file << Buffer << std::endl;

	FormatAGCUnitVector(Buffer2, PHASEB / PI2);
	sprintf_s(Buffer, 255, "       2DEC  %s      # BARG", Buffer2);
	file << Buffer << std::endl;

	FormatAGCUnitVector(Buffer2, OMEGAB / PI2);
	sprintf_s(Buffer, 255, "       2DEC* %s B+1* # 1/32", Buffer2);
	file << Buffer << std::endl;

	FormatAGCUnitVector(Buffer2, CMOD / PI2);
	sprintf_s(Buffer, 255, "       2DEC* %s B+1* # CMOD", Buffer2);
	file << Buffer << std::endl;

	FormatAGCUnitVector(Buffer2, CARG / PI2);
	sprintf_s(Buffer, 255, "       2DEC  %s      # CARG", Buffer2);
	file << Buffer << std::endl;

	FormatAGCUnitVector(Buffer2, OMEGAC / PI2);
	sprintf_s(Buffer, 255, "       2DEC* %s B+1* # 1/365", Buffer2);
	file << Buffer << std::endl;

	file << std::endl << "Solar Ephemeris error over time:" << std::endl;

	MJD = MJD_0;

	VECTOR3 U_ES1, U_ES2, R_ES, V_ES;
	double t;

	while (MJD < MJD_0 + 370.0)
	{
		//Time since July 1st, midnight
		t = MJD - MJD_0;
		//Sun unit vector with LGC ephemeris
		U_ES1 = EvaluateSolarEphemeris(t, LOSO, LOSR, CMOD, OMEGAC, CARG, K1, K3);
		//Sun unit vector with Orbiter ephemeris
		agcCelBody_RH(BODY_EARTH, MJD, EPHEM_TRUEPOS | EPHEM_TRUEVEL, &R_ES, &V_ES);
		R_ES = mul(Rot, -R_ES);
		U_ES2 = unit(R_ES);
		//Angle between vectors
		ang = acos(dotp(U_ES1, U_ES2));

		sprintf_s(Buffer, 256, "%.3f = %.8lf (degrees)", MJD, ang*DEG);
		file << Buffer << std::endl;

		MJD += 10.0;
	}

	file << std::endl << "Lunar Ephemeris error over time:" << std::endl;

	VECTOR3 U_EM1, U_EM2, R_EM;

	MJD = MJD_0;

	while (MJD < MJD_0 + 370.0)
	{
		//Time since July 1st, midnight
		t = MJD - MJD_0;
		//Sun unit vector with LGC ephemeris
		U_EM1 = EvaluateLunarEphemeris(t, LOMO, LOMR, A, OMEGAA, PHASEA, B, OMEGAB, PHASEB, LONO, LONR, K1, K2, K3, K4);
		//Moon unit vector with Orbiter ephemeris
		agcCelBody_RH(BODY_MOON, MJD, EPHEM_TRUEPOS, &R_EM);
		R_EM = mul(Rot, R_EM);
		U_EM2 = unit(R_EM);
		//Angle between vectors
		ang = acos(dotp(U_EM1, U_EM2));

		sprintf_s(Buffer, 256, "%.3f = %.8lf (degrees)", MJD, ang*DEG);
		file << Buffer << std::endl;

		MJD += 10.0;
	}

	file.close();
}