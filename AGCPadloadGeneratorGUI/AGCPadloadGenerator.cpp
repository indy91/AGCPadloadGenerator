#include "AGCPadloadGenerator.h"
#include "Vsop87.h"
#include <string>

void ELP82_init();
void ELP82_exit();
int  ELP82_read(double prec);
int  ELP82(double mjd, double *r);

#define SPS_THRUST					91188.544		// CMC fixed constant
#define SPS_ISP						 3080.0

const double R_Moon = 1.73809e6;			///< Radius of the moon
const double R_Earth = 6373338.0;			///< Radius of the Earth at the launch pad
const double LBS2KG = 0.45359237;			///< Pound mass to kilograms
const double LBF2N = LBS2KG * 9.80665;		///< Pound mass to kilograms
const double FT2M = 0.3048;					///< Feet to meters
const double w_Earth = 7.29211515e-05;
const double a_Fisher = 6373338.0;		//Semi-major axis of Fisher ellipsoid, should be 6378166.0
const double b_Fisher = 6373338.0;		//Semi-minor axis of Fisher ellipsoid, should be 6356784.0

void mjddate(double mjd, int &year, int &month, int &day)
{
	//Convert MJD to Year, Month, Day

	double ijd, c, e, h;
	int a, b, f;

	h = 24.0 * modf(mjd, &ijd);
	if (ijd < -100840) {
		c = ijd + 2401525.0;
	}
	else {
		b = (int)((ijd + 532784.75) / 36524.25);
		c = ijd + 2401526.0 + (b - b / 4);
	}
	a = (int)((c - 122.1) / 365.25);
	e = 365.0 * a + a / 4;
	f = (int)((c - e) / 30.6001);

	day = (int)(c - e + 0.5) - (int)(30.6001*f);
	month = f - 1 - 12 * (f / 14);
	year = a - 4715 - ((7 + month) / 10);
}

int SingleToBuffer(double x, double SF, bool TwosComplement = false)
{
	double x2;
	int c;

	x2 = fabs(x) * pow(2, -SF + 14);

	c = (unsigned)round(x2);

	if (TwosComplement == false && c > 037777)
	{
		c = 037777;
	}

	if (x < 0.0)
	{
		// Polarity change
		c = 0x7FFF & (~c);
		if (TwosComplement)
		{
			c++;
		}
	}

	return c;
}

void DoubleToBuffer(double x, int SF, int &c1, int &c2)
{
	double x2;

	x2 = fabs(x) * pow(2, -SF + 14);

	c1 = (unsigned)x2;
	x2 = x2 - (double)c1;
	x2 *= pow(2, 14);
	c2 = (unsigned)round(x2);
	if (c2 > 037777)
	{
		c2 = 037777;
	}

	if (x < 0.0) c1 = 0x7FFF & (~c1); // Polarity change
	if (x < 0.0) c2 = 0x7FFF & (~c2); // Polarity change
}

void TripleToBuffer(double x, int SF, int &c1, int &c2, int &c3)
{
	double x2;

	x2 = fabs(x) * pow(2, -SF + 14);

	c1 = (unsigned)x2;
	x2 = x2 - (double)c1;
	x2 *= pow(2, 14);
	c2 = (unsigned)x2;
	x2 = x2 - (double)c2;
	x2 *= pow(2, 14);
	c3 = (unsigned)round(x2);
	if (c3 > 037777)
	{
		c3 = 037777;
	}

	if (x < 0.0) c1 = 0x7FFF & (~c1); // Polarity change
	if (x < 0.0) c2 = 0x7FFF & (~c2); // Polarity change
	if (x < 0.0) c3 = 0x7FFF & (~c3); // Polarity change
}

void agcPolar2Cartesian_LH(VECTOR3 *_pv, VECTOR3 *_vv)
{
	double cx = cos((*_pv).x), sx = sin((*_pv).x);
	double cy = sin((*_pv).y), sy = cos((*_pv).y);
	double r = (*_pv).z * AU;

	VECTOR3 _r = _V(sy*cx, cy, sy*sx);

	*_pv = _r * r;

	if (_vv == NULL) return;

	VECTOR3 _y = _V(-cy * cx, sy, -cy * sx);
	VECTOR3 _x = _V(-sx, 0, cx);

	*_vv = _r * ((*_vv).z*AU) + _y * (r*(*_vv).y) + _x * (r*(*_vv).x*sy);
}

bool SolveSystem(int n, double *A, double *b, double *x, double *det)
{
	int e = 0, *p = new int[n];
	for (int i = 0;i < n;i++) p[i] = i;
	for (int k = 0;k < n;k++) {
		int r = 0; double d = 0.0;
		for (int s = k;s < n;s++) if (fabs(A[s*n + k]) > d) { d = fabs(A[s*n + k]); r = s; }
		if (d == 0.0) { delete[]p; return false; }
		if (r != k) { // Do Swaps
			for (int i = 0;i < n;i++) { double x = A[k*n + i]; A[k*n + i] = A[r*n + i]; A[r*n + i] = x; }
			int x = p[k]; p[k] = p[r]; p[r] = x; e++;
		}
		for (int i = k + 1;i < n;i++) { A[i*n + k] /= A[k*n + k]; for (int j = k + 1;j < n;j++) A[i*n + j] -= (A[i*n + k] * A[k*n + j]); }
	}
	for (int i = 0;i < n;i++) { x[i] = b[p[i]];	for (int j = 0;j < i;j++) x[i] -= A[i*n + j] * x[j]; }
	for (int i = n - 1;i >= 0;i--) { for (int j = i + 1;j < n;j++) x[i] -= A[i*n + j] * x[j]; x[i] /= A[i*n + i]; }
	if (det) { *det = 1.0; for (int i = 0;i < n;i++) *det *= A[i*n + i]; if (e & 1) *det *= -1.0; }
	delete[]p;
	return true;
}

bool SolveSeries(double *x, double *y, int ndata, double *out, int m)
{

	double *v = new double[m];
	double *q = new double[m];
	double *M = new double[m*m];

	memset(M, 0, m*m * sizeof(double));
	memset(q, 0, m * sizeof(double));

	for (int i = 0;i < ndata;i++) {
		v[0] = 1.0;
		for (int f = 1;f < m;f++) v[f] = v[f - 1] * x[i];
		for (int f = 0;f < m;f++) for (int g = 0;g < m;g++) M[f*m + g] += (v[f] * v[g]);
		for (int f = 0;f < m;f++) q[f] += (v[f] * y[i]);
	}

	bool bRet = SolveSystem(m, M, q, out, NULL);
	delete[]v; delete[]q; delete[]M;
	return bRet;
}

IMUBiasCompensationData::IMUBiasCompensationData()
{
	PBIASX = 0.0;
	PIPASCFX = 0.0;
	PBIASY = 0.0;
	PIPASCFY = 0.0;
	PBIASZ = 0.0;
	PIPASCFZ = 0.0;
	NBDX = 0.0;
	NBDY = 0.0;
	NBDZ = 0.0;
	ADIAX = 0.0;
	ADIAY = 0.0;
	ADIAZ = 0.0;
	ADSRAX = 0.0;
	ADSRAY = 0.0;
	ADSRAZ = 0.0;
}

AGCPadloadGenerator::AGCPadloadGenerator()
{
	ELP82_init();
	earth.clbkInit();

	LaunchAzimuth = 72.0;
	RPVAR = 4000000.0;
	S22WSUBL = 10000.0;
	WORBPOS = 0.0;
	WORBVEL = 0.0;
	CDUCHKWD = 0.05;			// 5 cs
	DVTHRESH = 2.0*0.3048;		// 2 ft/s
	HORIZALT = 28000.0;			// 28000 m
	ALTVAR = 1.52168e-5;		// rad^2

	BLOCKI.T_ATL = 1400.0;
	BLOCKI.lat_ATL = 28.29028886;
	BLOCKI.lng_ATL = -19.5;
	BLOCKI.T_PAC = 30921.42;
	BLOCKI.lat_PAC = 30.04649677;
	BLOCKI.lng_PAC = -171.0;
	BLOCKI.e_SPS1 = 0.5934490037;
	BLOCKI.a_SPS1 = 15487553.0;
	BLOCKI.e_SPS2 = 0.999071629063702; // 0.999071629;
	BLOCKI.a_SPS2 = 6891085630.5;
	BLOCKI.TROLL = 9.0;
	BLOCKI.TPITCH = 10.7;
	BLOCKI.TENDPITCH = 134.0;

	BLOCKII.CSMMass = 0.0;
	BLOCKII.LMMass = 0.0;
	BLOCKII.TotalMass = 0.0;
	BLOCKII.WRENDPOS = 10000.0;	// 10000 ft
	BLOCKII.WRENDVEL = 10.0;	// 10 ft/s
	BLOCKII.WSHAFT = 15.0; //mr
	BLOCKII.WTRUN = 15.0; //mr
	BLOCKII.RMAX = 2000.0;		// 2000 ft
	BLOCKII.VMAX = 2.0;			// 2 ft/s
	BLOCKII.RVARMIN = 40000.0; // (200 ft)^2
	BLOCKII.SHAFTVAR = 1.0;
	BLOCKII.TRUNVAR = 1.0;
	BLOCKII.WSURFPOS = 0.0;
	BLOCKII.WSURFVEL = 0.0;
	BLOCKII.HIASCENT = 10900.0; //lbs
	BLOCKII.AGSK = 100.0; //hours
	BLOCKII.ROLLTIME = 6.0; //deg
	BLOCKII.PITCHTIME = 6.0; //deg
	BLOCKII.TLAND = 0.0;
	BLOCKII.LAT_SPL = 26.48; //Target for 72° azimuth
	BLOCKII.LNG_SPL = -17.05;
	BLOCKII.PACTOFF = 0.0;
	BLOCKII.YACTOFF = 0.0;
	BLOCKII.LADPAD = 0.3;
	BLOCKII.LODPAD = 0.18;
	BLOCKII.ALFAPAD = -20.0;
	BLOCKII.P37RANGE = 0.0;
	for (int i = 0; i < 7; i++)
	{
		BLOCKII.POLYNUM[i] = 0.0;
	}
	BLOCKII.RPSTART = 0.0;
	BLOCKII.POLYSTOP = 0.0;

	for (int i = 0; i < 5; i++)
	{
		BLOCKII.ABSC[i] = BLOCKII.SLOPE[i] = 0.0;
	}
	BLOCKII.J1PARM = BLOCKII.K1PARM = BLOCKII.J2PARM = BLOCKII.K2PARM = BLOCKII.THETCRIT = BLOCKII.RAMIN = 0.0;
	BLOCKII.DELQFIX = 100.0;
	BLOCKII.IGNAOSQ = 0.0;
	BLOCKII.IGNAOSR = 0.0;
	BLOCKII.DELTTFAP = -90.0;

	CMCDATA.EIMP1SEC = 19575.168;
	CMCDATA.EFIMP01 = 24484.268;
	CMCDATA.EFIMP16 = 20281.0;
	CMCDATA.WMIDPOS = 30000.0;	// 10000 ft
	CMCDATA.WMIDVEL = 30.0;		// 10 ft/s
}

AGCPadloadGenerator::~AGCPadloadGenerator()
{
	ELP82_exit();
}

void AGCPadloadGenerator::SaveEMEM(int address, int value)
{
	EMEM temp;

	temp.address = address;
	temp.value = value;

	arr.push_back(temp);
}

void AGCPadloadGenerator::WriteEMEM(int address, int value, bool cmc)
{
	if (cmc)
	{
		sprintf_s(Buffer, "  CMPAD%04o %o", address, value);
	}
	else
	{
		sprintf_s(Buffer, "  LMPAD %04o %05o", address, value);
	}
	myfile << Buffer << std::endl;
}

int AGCPadloadGenerator::RunLGC()
{
	//Clear erasable data
	arr.clear();
	//Open padload file
	myfile.open("Padload.txt");

	//Launch MJD rounded down to the next 12 hours
	double AGCEphemStartTime = floor(LaunchMJD*2.0) / 2.0;
	//Midnight day before launch
	PrelaunchMJD = floor(LaunchMJD) - 1.0;

	//Get rope number
	AGCVersions AGCVersion = GetLGCVersion(RopeName);
	if (AGCVersion == AGCVersions::AGCVersionError) return 1;

	//Get PIOS data set

	if (PIOSDataSetName == "") return 1;

	PIOSDataSet DataSet;
	int err = GetPIOSDataSet(PIOSDataSetName, DataSet);
	if (err) return 4;

	switch (AGCVersion)
	{
		case Sundance306:
			Sundance306Defaults();
			break;
		case Luminary069:
			Luminary069Padload(false);
			break;
		case Luminary069R2:
			Luminary069Padload(true);
			break;
		case Luminary099:
			Luminary099Padload();
			break;
		case Luminary116:
			Luminary116Padload();
			break;
		case Luminary131:
			return 2; //Error
		case Luminary131R1:
			Luminary131Padload();
			break;
		case Luminary178:
			Luminary178Padload();
			break;
		case Luminary210:
			Luminary210Padload();
			break;
		case Zerlina56:
			Zerlina56Padload();
			break;
		default: //Error
			return 2;
	}

	//Check on TEPHEM (possible TEPHEM limit or too early)
	double limit = DataSet.ExtendedLimit ? 963.0 : 481.0;
	if (LaunchMJD - DataSet.t0 > limit || LaunchMJD < DataSet.t0)
	{
		return 3;
	}

	//Calculate padload TEPHEM
	TEPHEM = (LaunchMJD - DataSet.t0)*24.0*3600.0*100.0;

	TripleToBuffer(TEPHEM, 42, iTemp, iTemp2, iTemp3);
	SaveEMEM(01706, iTemp);
	SaveEMEM(01707, iTemp2);
	SaveEMEM(01710, iTemp3);

	AGCCorrectionVectors(DataSet, AGCEphemStartTime, T_UNITW, T_504LM, false, RopeName == "Sundance306");

	//End
	std::sort(arr.begin(), arr.end());

	//Count
	sprintf_s(Buffer, "  LMPADCNT %d", arr.size());
	myfile << Buffer << std::endl;

	for (unsigned i = 0;i < arr.size();i++)
	{
		WriteEMEM(arr[i].address, arr[i].value, false);
	}

	myfile.close();

	return 0;
}

void AGCPadloadGenerator::RunBlockI()
{
	arr.clear();

	if (Pad == "LC-39A")
	{
		SetPadData(Launchpad::LC39A);
	}
	else if (Pad == "LC-34")
	{
		SetPadData(Launchpad::LC34);
	}
	else
	{
		return;
	}

	//MJD of preceding July 1st, midnight
	double AGCEphemTEphemZero, A_Z0;
	int Epoch;

	if (RopeName == "Corona261")
	{
		AGCEphemTEphemZero = 39307.0;
		Epoch = 1967;
		A_Z0 = HANGLE(Epoch, Epoch - 1, 182); //July 1st
	}
	else if (RopeName == "Solarium055")
	{
		AGCEphemTEphemZero = 39672.0;
		Epoch = 1968;
		A_Z0 = HANGLE(Epoch, Epoch - 1, 182); //July 1st. 278.368527*RAD;
	}
	else
	{
		return;
	}

	double brcsmjd = MJDOfNBYEpoch(Epoch);
	PrelaunchMJD = floor(LaunchMJD);

	BLOCKI.DTEPOCH = Solarium055DTEPOCHCalculation(A_Z0, AGCEphemTEphemZero, PrelaunchMJD, PadLong*RAD);

	if (RopeName == "Corona261")
	{
		CoronaDefaults();
	}
	else
	{
		SolariumDefaults();
	}

	myfile.open("Padload.txt");

	//End
	std::sort(arr.begin(), arr.end());

	for (unsigned i = 0; i < arr.size(); i++)
	{
		WriteEMEM(arr[i].address, arr[i].value, true);
	}

	myfile.close();
}

int AGCPadloadGenerator::RunCMC()
{
	arr.clear();

	if (Pad == "LC-39A")
	{
		SetPadData(Launchpad::LC39A);
	}
	else if (Pad == "LC-34")
	{
		SetPadData(Launchpad::LC34);
	}
	else
	{
		SetPadData(Launchpad::LC39B);
	}

	myfile.open("Padload.txt");

	//Get rope number
	AGCVersions AGCVersion = GetCMCVersion(RopeName);
	if (AGCVersion == AGCVersions::AGCVersionError) return 1;

	//Get PIOS data set
	if (PIOSDataSetName == "") return 1;

	PIOSDataSet DataSet;
	int err = GetPIOSDataSet(PIOSDataSetName, DataSet);
	if (err) return 4;

	//Launch MJD rounded down to the next 12 hours
	double AGCEphemStartTime = floor(LaunchMJD*2.0) / 2.0;
	//Midnight day before launch
	PrelaunchMJD = floor(LaunchMJD) - 1.0;

	switch (AGCVersion)
	{
	case Colossus237:
		Colossus237_249_Defaults(false);
		break;
	case Colossus249:
		Colossus237_249_Defaults(true);
		break;
	case Comanche045:
		Comanche45Padload(true);
		break;
	case Comanche055:
		Comanche55Padload();
		break;
	case Comanche067:
		Comanche67Padload();
		break;
	case Comanche072:
		Comanche72Padload();
		break;
	case Comanche108:
		Comanche108Padload();
		break;
	case Artemis072:
		Artemis72Padload();
		break;
	case Skylark048:
	{
		int Year, Month, Day;

		//Get Year and Month of prelaunch
		mjddate(PrelaunchMJD, Year, Month, Day);

		//Get Year of July 1st that preceeds the PrelaunchMJD
		if (Month < 7)
		{
			Year--;
		}

		//Get MJD at midnight July 1st that preceeds PrelaunchMJD
		DataSet.t0 = JD2MJD(TJUDAT(Year, 7, 1));

		Skylark048Padload();
		SkylarkSolarEphemeris(MJD2JD(LaunchMJD), MJD2JD(DataSet.t0));
		SkylarkCorrectionMatrix(MJD2JD(LaunchMJD), MJD2JD(DataSet.t0));
	}
	break;
	default:
		return 2; //Error
	}

	//Check on TEPHEM (possible TEPHEM limit or too early)
	double limit = DataSet.ExtendedLimit ? 963.0 : 481.0;
	if (LaunchMJD - DataSet.t0 > limit || LaunchMJD < DataSet.t0)
	{
		return 3;
	}

	//Calculate padload TEPHEM
	TEPHEM = (PrelaunchMJD - DataSet.t0)*24.0*3600.0*100.0;

	TripleToBuffer(TEPHEM, 42, iTemp, iTemp2, iTemp3);

	if (AGCVersion == Skylark048)
	{
		SaveEMEM(01700, iTemp);
		SaveEMEM(01701, iTemp2);
		SaveEMEM(01702, iTemp3);
	}
	else
	{
		SaveEMEM(01706, iTemp);
		SaveEMEM(01707, iTemp2);
		SaveEMEM(01710, iTemp3);

		AGCCorrectionVectors(DataSet, AGCEphemStartTime, T_UNITW, T_504LM, true, false);
		AGCEphemeris(AGCEphemStartTime, DataSet.epoch, DataSet.t0, EphemerisSpan);
	}

	//End
	std::sort(arr.begin(), arr.end());

	for (unsigned i = 0;i < arr.size();i++)
	{
		WriteEMEM(arr[i].address, arr[i].value, true);
	}

	myfile.close();

	return 0;
}

void AGCPadloadGenerator::SetPadData(Launchpad pad)
{
	if (pad == Launchpad::LC39A)
	{
		PadLat = 28.60842218;
		PadLong = 279.3958666;
		PadAlt = 89.4;
		PadAzimuth = -90.0;

		TAZEL[0] = 253.192 - 360.0;
		TAZEL[1] = -2.0112;
		TAZEL[2] = -68.81059444;// 291.187 - 360.0;
		TAZEL[3] = -2.0158;
	}
	else if (pad == Launchpad::LC39B)
	{
		PadLat = 28.62687861;
		PadLong = 279.3789091;
		PadAlt = 89.4;
		PadAzimuth = -90.0;

		TAZEL[0] = -106.7757496;
		TAZEL[1] = -1.757722222;
		TAZEL[2] = -63.7823055;
		TAZEL[3] = -1.678888922;
	}
	else if (pad == Launchpad::LC34)
	{
		PadLat = 28.5217969;
		PadLong = 279.4387535;
		PadAlt = 44.0;
		PadAzimuth = -80.0;

		TAZEL[0] = 309.6 - 360.0;
		TAZEL[1] = -0.1;
		TAZEL[2] = 270.0 - 360.0;
		TAZEL[3] = 0.0;
	}
	else
	{

	}
}

void AGCPadloadGenerator::LGCDefaults(bool EarlyPIPABias, bool mass)
{
	//Contains padloads that never change their address in all of Sundance and Luminary

	IMUCompensation(false, EarlyPIPABias);

	//GCOMPSW
	SaveEMEM(01477, 0);

	//X789
	SaveEMEM(01700, 0);
	SaveEMEM(01701, 0);
	SaveEMEM(01702, 0);
	SaveEMEM(01703, 0);
	SaveEMEM(01704, 0);
	SaveEMEM(01705, 0);

	//WRENDPOS
	iTemp = SingleToBuffer(BLOCKII.WRENDPOS*0.3048, 14);
	SaveEMEM(02000, iTemp);
	//WRENDVEL
	iTemp = SingleToBuffer(BLOCKII.WRENDVEL*0.3048 / 100.0, 0);
	SaveEMEM(02001, iTemp);
	//WSHAFT
	iTemp = SingleToBuffer(BLOCKII.WSHAFT / 1000.0, -5);
	SaveEMEM(02002, iTemp);
	//WTRUN
	iTemp = SingleToBuffer(BLOCKII.WTRUN / 1000.0, -5);
	SaveEMEM(02003, iTemp);
	//RMAX
	iTemp = SingleToBuffer(BLOCKII.RMAX*0.3048, 19);
	SaveEMEM(02004, iTemp);
	//VMAX
	iTemp = SingleToBuffer(BLOCKII.VMAX*0.3048 / 100.0, 7);
	SaveEMEM(02005, iTemp);

	//HIASCENT
	iTemp = SingleToBuffer(BLOCKII.HIASCENT*LBS2KG, 16);
	SaveEMEM(03000, iTemp);

	//ROLLTIME
	iTemp = SingleToBuffer(BLOCKII.ROLLTIME / 0.2*100.0, 14); //6° converted to 0.2°/s speed
	SaveEMEM(03001, iTemp);
	//PITTIME
	iTemp = SingleToBuffer(BLOCKII.PITCHTIME / 0.2*100.0, 14); //6° converted to 0.2°/s speed
	SaveEMEM(03002, iTemp);
	//DKTRAP
	iTemp = SingleToBuffer(-1.4 / 360.0, -3);
	SaveEMEM(03003, iTemp);
	//DKOMEGAN
	iTemp = SingleToBuffer(10.0, 14);
	SaveEMEM(03004, iTemp);
	//DKKAOSN
	iTemp = SingleToBuffer(60.0, 14);
	SaveEMEM(03005, iTemp);
	//LMTRAP
	iTemp = SingleToBuffer(-1.4 / 360.0, -3);
	SaveEMEM(03006, iTemp);
	//LMOMEGAN
	iTemp = SingleToBuffer(0.0, 14);
	SaveEMEM(03007, iTemp);
	//LMKAOSN
	iTemp = SingleToBuffer(60.0, 14);
	SaveEMEM(03010, iTemp);
	//DKDB
	SaveEMEM(03011, 0200);

	//ATIGINC
	dTemp = 180.0*100.0;
	DoubleToBuffer(dTemp, 28, iTemp, iTemp2);
	SaveEMEM(03400, iTemp);
	SaveEMEM(03401, iTemp2);
	//PTIGINC
	DoubleToBuffer(dTemp, 28, iTemp, iTemp2);
	SaveEMEM(03402, iTemp);
	SaveEMEM(03403, iTemp2);
}

void AGCPadloadGenerator::Sundance306Defaults()
{
	LGCDefaults(true);

	//DUMPCNT
	SaveEMEM(0333, 010000);

	//PHSPRDT2 - restart protection during P00
	SaveEMEM(01057, 013000);

	//MASS
	DoubleToBuffer(BLOCKII.TotalMass*LBS2KG, 16, iTemp, iTemp2);
	SaveEMEM(01243, iTemp);
	SaveEMEM(01244, iTemp2);

	//LEMMASS
	iTemp = SingleToBuffer(BLOCKII.LMMass*LBS2KG, 16);
	SaveEMEM(01335, iTemp);
	//CSMMASS
	iTemp = SingleToBuffer(BLOCKII.CSMMass*LBS2KG, 16);
	SaveEMEM(01336, iTemp);

	//SHAFTVAR
	iTemp = SingleToBuffer(BLOCKII.SHAFTVAR*1e-6, -12);
	SaveEMEM(02006, iTemp);
	//TRUNVAR
	iTemp = SingleToBuffer(BLOCKII.TRUNVAR*1e-6, -12);
	SaveEMEM(02007, iTemp);

	//AGSK
	DoubleToBuffer(BLOCKII.AGSK*3600.0*100.0, 28, iTemp, iTemp2);
	SaveEMEM(02016, iTemp);
	SaveEMEM(02017, iTemp2);

	//RANGEVAR
	DoubleToBuffer(1.111111111e-5, -12, iTemp, iTemp2);
	SaveEMEM(02364, iTemp);
	SaveEMEM(02365, iTemp2);
	//RATEVAR
	DoubleToBuffer(1.877777000e-5, -12, iTemp, iTemp2);
	SaveEMEM(02366, iTemp);
	SaveEMEM(02367, iTemp2);
	//RVARMIN
	SaveEMEM(02370, 0);
	SaveEMEM(02371, 0);
	iTemp = SingleToBuffer(66.0, 12);
	SaveEMEM(02372, iTemp);
	//VVARMIN
	SaveEMEM(02373, 0);
	iTemp = SingleToBuffer(1.7445e-6, -12);
	SaveEMEM(02374, iTemp);

	//AOTAZ
	iTemp = SingleToBuffer(-60.0 / 360.0, -1, true);
	SaveEMEM(03404, iTemp);
	iTemp = SingleToBuffer(0.0 / 360.0, -1, true);
	SaveEMEM(03405, iTemp);
	iTemp = SingleToBuffer(60.0 / 360.0, -1, true);
	SaveEMEM(03406, iTemp);
	//AOTEL
	iTemp = SingleToBuffer(45.0 / 360.0, -1);
	SaveEMEM(03407, iTemp);
	iTemp = SingleToBuffer(45.0 / 360.0, -1);
	SaveEMEM(03410, iTemp);
	iTemp = SingleToBuffer(45.0 / 360.0, -1);
	SaveEMEM(03411, iTemp);

	//ZOOMTIME
	iTemp = SingleToBuffer(26.0*100.0, 14);
	SaveEMEM(03421, iTemp);
}

void AGCPadloadGenerator::Luminary069Padload(bool IsR2)
{
	LGCDefaults(true, true);

	//FLAGWRD3
	SaveEMEM(077, 02000);
	//FLAGWRD8
	SaveEMEM(0104, 0);
	//FLAGWRD10
	SaveEMEM(0106, 0);

	//DUMPCNT
	SaveEMEM(0333, 010000); //2 dumps

	//MASS
	DoubleToBuffer(BLOCKII.TotalMass*LBS2KG, 16, iTemp, iTemp2);
	SaveEMEM(01244, iTemp);
	SaveEMEM(01245, iTemp2);

	//LEMMASS
	iTemp = SingleToBuffer(BLOCKII.LMMass*LBS2KG, 16);
	SaveEMEM(01331, iTemp);
	//CSMMASS
	iTemp = SingleToBuffer(BLOCKII.CSMMass*LBS2KG, 16);
	SaveEMEM(01332, iTemp);

	//TETCSM
	SaveEMEM(01570, 037777);
	//TETLEM
	SaveEMEM(01642, 037777);

	//RANGEVAR
	DoubleToBuffer(1.111111111e-5, -12, iTemp, iTemp2);
	SaveEMEM(01770, iTemp);
	SaveEMEM(01771, iTemp2);
	//RATEVAR
	DoubleToBuffer(1.877777000e-5, -12, iTemp, iTemp2);
	SaveEMEM(01772, iTemp);
	SaveEMEM(01773, iTemp2);
	//RVARMIN
	iTemp = SingleToBuffer(66.0, 12);
	SaveEMEM(01774, iTemp);
	//VVARMIN
	iTemp = SingleToBuffer(1.7445e-6, -12);
	SaveEMEM(01775, iTemp);

	//WSURFPOS
	iTemp = SingleToBuffer(BLOCKII.WSURFPOS*FT2M, 14);
	SaveEMEM(02006, iTemp);
	//WSURFVEL
	iTemp = SingleToBuffer(BLOCKII.WSURFVEL*FT2M / 100.0, 0);
	SaveEMEM(02007, iTemp);

	//SHAFTVAR
	iTemp = SingleToBuffer(BLOCKII.SHAFTVAR*1e-6, -12);
	SaveEMEM(02010, iTemp);
	//TRUNVAR
	iTemp = SingleToBuffer(BLOCKII.TRUNVAR*1e-6, -12);
	SaveEMEM(02011, iTemp);

	//AGSK
	DoubleToBuffer(BLOCKII.AGSK*3600.0*100.0, 28, iTemp, iTemp2);
	SaveEMEM(02020, iTemp);
	SaveEMEM(02021, iTemp2);

	//RLS
	RLS = r_from_latlong(LSLat*RAD, LSLng*RAD, LSAlt, BODY_MOON, 1);
	DoubleToBuffer(RLS.x, 27, iTemp, iTemp2);
	SaveEMEM(02022, iTemp);
	SaveEMEM(02023, iTemp2);
	DoubleToBuffer(RLS.y, 27, iTemp, iTemp2);
	SaveEMEM(02024, iTemp);
	SaveEMEM(02025, iTemp2);
	DoubleToBuffer(RLS.z, 27, iTemp, iTemp2);
	SaveEMEM(02026, iTemp);
	SaveEMEM(02027, iTemp2);

	//For now hardcoded
	BLOCKII.RBRFG.x = 171.835;
	BLOCKII.RBRFG.y = 0.0;
	BLOCKII.RBRFG.z = -10678.596;

	BLOCKII.VBRFG.x = -105.876;
	BLOCKII.VBRFG.y = 0.0;
	BLOCKII.VBRFG.z = -1.04;

	BLOCKII.ABRFG.x = 0.6241;
	BLOCKII.ABRFG.y = 0.0;
	BLOCKII.ABRFG.z = -9.1044;

	BLOCKII.VBRFG_star = -2.34;
	BLOCKII.ABRFG_star = -54.6264;
	BLOCKII.JBRFG_star = -0.15061416;

	BLOCKII.RAPFG.x = 111.085;
	BLOCKII.RAPFG.y = 0.0;
	BLOCKII.RAPFG.z = -26.794;

	BLOCKII.VAPFG.x = -4.993;
	BLOCKII.VAPFG.y = 0.0;
	BLOCKII.VAPFG.z = 0.248;

	BLOCKII.AAPFG.x = -0.2624;
	BLOCKII.AAPFG.y = 0.0;
	BLOCKII.AAPFG.z = -0.512;

	BLOCKII.VAPFG_star = 0.558;
	BLOCKII.AAPFG_star = -3.072;
	BLOCKII.JAPFG_star = 0.01446176;

	//TLAND
	DoubleToBuffer(BLOCKII.TLAND*3600.0*100.0, 28, iTemp, iTemp2);
	SaveEMEM(02400, iTemp);
	SaveEMEM(02401, iTemp2);

	//TARGRDG
	DoubleToBuffer(BLOCKII.RBRFG.x*FT2M, 24, iTemp, iTemp2);
	SaveEMEM(02402, iTemp);
	SaveEMEM(02403, iTemp2);
	DoubleToBuffer(BLOCKII.RBRFG.y*FT2M, 24, iTemp, iTemp2);
	SaveEMEM(02404, iTemp);
	SaveEMEM(02405, iTemp2);
	DoubleToBuffer(BLOCKII.RBRFG.z*FT2M, 24, iTemp, iTemp2);
	SaveEMEM(02406, iTemp);
	SaveEMEM(02407, iTemp2);

	//TARGVDG
	DoubleToBuffer(BLOCKII.VBRFG.x*FT2M / 100.0, 10, iTemp, iTemp2);
	SaveEMEM(02410, iTemp);
	SaveEMEM(02411, iTemp2);
	DoubleToBuffer(BLOCKII.VBRFG.y*FT2M / 100.0, 10, iTemp, iTemp2);
	SaveEMEM(02412, iTemp);
	SaveEMEM(02413, iTemp2);
	DoubleToBuffer(BLOCKII.VBRFG.z*FT2M / 100.0, 10, iTemp, iTemp2);
	SaveEMEM(02414, iTemp);
	SaveEMEM(02415, iTemp2);

	//TARGADG
	DoubleToBuffer(BLOCKII.ABRFG.x*FT2M / 10000.0, -4, iTemp, iTemp2);
	SaveEMEM(02416, iTemp);
	SaveEMEM(02417, iTemp2);
	DoubleToBuffer(BLOCKII.ABRFG.y*FT2M / 10000.0, -4, iTemp, iTemp2);
	SaveEMEM(02420, iTemp);
	SaveEMEM(02421, iTemp2);
	DoubleToBuffer(BLOCKII.ABRFG.z*FT2M / 10000.0, -4, iTemp, iTemp2);
	SaveEMEM(02422, iTemp);
	SaveEMEM(02423, iTemp2);

	//TTFVDGV
	DoubleToBuffer(BLOCKII.VBRFG_star*FT2M / 100.0, 13, iTemp, iTemp2);
	SaveEMEM(02424, iTemp);
	SaveEMEM(02425, iTemp2);
	//TTFADGZ
	DoubleToBuffer(BLOCKII.ABRFG_star*FT2M / 10000.0, -4, iTemp, iTemp2);
	SaveEMEM(02426, iTemp);
	SaveEMEM(02427, iTemp2);
	//TTFJDGZ
	DoubleToBuffer(BLOCKII.JBRFG_star*FT2M / 1000000.0, -21, iTemp, iTemp2);
	SaveEMEM(02430, iTemp);
	SaveEMEM(02431, iTemp2);

	//TARGRDG24
	DoubleToBuffer(BLOCKII.RAPFG.x*FT2M, 24, iTemp, iTemp2);
	SaveEMEM(02432, iTemp);
	SaveEMEM(02433, iTemp2);
	DoubleToBuffer(BLOCKII.RAPFG.y*FT2M, 24, iTemp, iTemp2);
	SaveEMEM(02434, iTemp);
	SaveEMEM(02435, iTemp2);
	DoubleToBuffer(BLOCKII.RAPFG.z*FT2M, 24, iTemp, iTemp2);
	SaveEMEM(02436, iTemp);
	SaveEMEM(02437, iTemp2);

	//TARGVDG24
	DoubleToBuffer(BLOCKII.VAPFG.x*FT2M / 100.0, 10, iTemp, iTemp2);
	SaveEMEM(02440, iTemp);
	SaveEMEM(02441, iTemp2);
	DoubleToBuffer(BLOCKII.VAPFG.y*FT2M / 100.0, 10, iTemp, iTemp2);
	SaveEMEM(02442, iTemp);
	SaveEMEM(02443, iTemp2);
	DoubleToBuffer(BLOCKII.VAPFG.z*FT2M / 100.0, 10, iTemp, iTemp2);
	SaveEMEM(02444, iTemp);
	SaveEMEM(02445, iTemp2);

	//TARGADG24
	DoubleToBuffer(BLOCKII.AAPFG.x*FT2M / 10000.0, -4, iTemp, iTemp2);
	SaveEMEM(02446, iTemp);
	SaveEMEM(02447, iTemp2);
	DoubleToBuffer(BLOCKII.AAPFG.y*FT2M / 10000.0, -4, iTemp, iTemp2);
	SaveEMEM(02450, iTemp);
	SaveEMEM(02451, iTemp2);
	DoubleToBuffer(BLOCKII.AAPFG.z*FT2M / 10000.0, -4, iTemp, iTemp2);
	SaveEMEM(02452, iTemp);
	SaveEMEM(02453, iTemp2);

	//TTFVDGZ24
	DoubleToBuffer(BLOCKII.VAPFG_star*FT2M / 100.0, 13, iTemp, iTemp2);
	SaveEMEM(02454, iTemp);
	SaveEMEM(02455, iTemp2);
	//TTFADGZ24
	DoubleToBuffer(BLOCKII.AAPFG_star*FT2M / 10000.0, -4, iTemp, iTemp2);
	SaveEMEM(02456, iTemp);
	SaveEMEM(02457, iTemp2);
	//TTFJDGZ24
	DoubleToBuffer(BLOCKII.JAPFG_star*FT2M / 1000000.0, -21, iTemp, iTemp2);
	SaveEMEM(02460, iTemp);
	SaveEMEM(02461, iTemp2);

	//DESIGNV
	DoubleToBuffer(BLOCKII.VIGN*FT2M / 100.0, 10, iTemp, iTemp2);
	SaveEMEM(02462, iTemp);
	SaveEMEM(02463, iTemp2);
	//DESIGNRX
	DoubleToBuffer(BLOCKII.RIGNX*0.3048, 24, iTemp, iTemp2);
	SaveEMEM(02464, iTemp);
	SaveEMEM(02465, iTemp2);
	//DESIGNRZ
	DoubleToBuffer(BLOCKII.RIGNZ*0.3048, 24, iTemp, iTemp2);
	SaveEMEM(02466, iTemp);
	SaveEMEM(02467, iTemp2);
	//DESKIGNX
	DoubleToBuffer(BLOCKII.KIGNXB4, 4, iTemp, iTemp2);
	SaveEMEM(02470, iTemp);
	SaveEMEM(02471, iTemp2);
	//DESKIGNY
	DoubleToBuffer(BLOCKII.KIGNYB8 / 0.3048, -16, iTemp, iTemp2);
	SaveEMEM(02472, iTemp);
	SaveEMEM(02473, iTemp2);
	//DESKIGNV
	DoubleToBuffer(BLOCKII.KIGNVB4*100.0, 18, iTemp, iTemp2);
	SaveEMEM(02474, iTemp);
	SaveEMEM(02475, iTemp2);
	//LOWCRIT
	SaveEMEM(02476, 04245);
	//HIGHCRIT
	SaveEMEM(02477, 04616);

	//DELQFIX
	DoubleToBuffer(200.0*0.3048, 24, iTemp, iTemp2);
	SaveEMEM(02500, iTemp);
	SaveEMEM(02501, iTemp2);

	//TBRKPNT
	iTemp = SingleToBuffer(600.0 *100.0, 17);
	SaveEMEM(02502, iTemp);

	//ABTVINJ1
	DoubleToBuffer(5605.0*FT2M / 100.0, 7, iTemp, iTemp2);
	SaveEMEM(02503, iTemp);
	SaveEMEM(02504, iTemp2);

	//ABTVINJ2
	DoubleToBuffer(5510.0*FT2M / 100.0, 7, iTemp, iTemp2);
	SaveEMEM(02505, iTemp);
	SaveEMEM(02506, iTemp2);

	//DOWNTORK
	SaveEMEM(03113, 0);
	SaveEMEM(03114, 0);
	SaveEMEM(03115, 0);
	SaveEMEM(03116, 0);
	SaveEMEM(03117, 0);
	SaveEMEM(03120, 0);

	//AOTAZ
	iTemp = SingleToBuffer(-60.0 / 360.0, -1, true);
	SaveEMEM(03404, iTemp);
	iTemp = SingleToBuffer(0.0 / 360.0, -1, true);
	SaveEMEM(03405, iTemp);
	iTemp = SingleToBuffer(60.0 / 360.0, -1, true);
	SaveEMEM(03406, iTemp);
	//AOTEL
	iTemp = SingleToBuffer(45.0 / 360.0, -1);
	SaveEMEM(03407, iTemp);
	iTemp = SingleToBuffer(45.0 / 360.0, -1);
	SaveEMEM(03410, iTemp);
	iTemp = SingleToBuffer(45.0 / 360.0, -1);
	SaveEMEM(03411, iTemp);

	//LRALPHA
	iTemp = SingleToBuffer(6.0 / 360.0, -1);
	SaveEMEM(03412, iTemp);
	//LRBETA1
	iTemp = SingleToBuffer(24.0 / 360.0, -1);
	SaveEMEM(03413, iTemp);
	//LRALPHA2
	iTemp = SingleToBuffer(6.0 / 360.0, -1);
	SaveEMEM(03414, iTemp);
	//LRBETA2
	iTemp = SingleToBuffer(0.0 / 360.0, -1);
	SaveEMEM(03415, iTemp);
	//LRHMAX
	iTemp = SingleToBuffer(50000.0*0.3048, 14);
	SaveEMEM(03416, iTemp);
	//LRVMAX
	iTemp = SingleToBuffer(2000.0*0.3048 / 100.0, 7);
	SaveEMEM(03417, iTemp);
	//LRWH
	iTemp = SingleToBuffer(0.35, -1);
	SaveEMEM(03420, iTemp);
	//LRWVZ
	iTemp = SingleToBuffer(0.3, 0);
	SaveEMEM(03421, iTemp);
	//LRWVY
	iTemp = SingleToBuffer(0.3, 0);
	SaveEMEM(03422, iTemp);
	//LRWVX
	iTemp = SingleToBuffer(0.3, 0);
	SaveEMEM(03423, iTemp);
	//ZOOMTIME
	iTemp = SingleToBuffer(26.0*100.0, 14);
	SaveEMEM(03424, iTemp);
	//TENDBRAK
	iTemp = SingleToBuffer(62.0*100.0, 17);
	SaveEMEM(03425, iTemp);
	//TENDAPPR
	iTemp = SingleToBuffer(12.0*100.0, 17);
	SaveEMEM(03426, iTemp);
	//RPCRTIME
	iTemp = SingleToBuffer(62.0*100.0, 17);
	SaveEMEM(03427, iTemp);
	//RPCRTQSW
	iTemp = SingleToBuffer(-0.5, 1);
	SaveEMEM(03430, iTemp);
}

void AGCPadloadGenerator::Luminary099Padload()
{
	LGCDefaults(true, true);
	Luminary099_116_Defaults();

	//FLAGWRD3
	SaveEMEM(077, 02000);

	//DUMPCNT
	SaveEMEM(0333, 010000); //2 dumps

	//MASS
	DoubleToBuffer(BLOCKII.TotalMass*LBS2KG, 16, iTemp, iTemp2);
	SaveEMEM(01244, iTemp);
	SaveEMEM(01245, iTemp2);

	//LEMMASS
	iTemp = SingleToBuffer(BLOCKII.LMMass*LBS2KG, 16);
	SaveEMEM(01331, iTemp);
	//CSMMASS
	iTemp = SingleToBuffer(BLOCKII.CSMMass*LBS2KG, 16);
	SaveEMEM(01332, iTemp);

	R2Model(01352);

	//RADSKAL
	SaveEMEM(01354, 0);
	SaveEMEM(01355, 0);
	//SKALSKAL
	SaveEMEM(01356, 0);

	//REFSMMAT
	SaveEMEM(01733, 012704);
	SaveEMEM(01734, 06264);
	SaveEMEM(01735, 012562);
	SaveEMEM(01736, 010723);
	SaveEMEM(01737, 01112);
	SaveEMEM(01740, 025001);

	//For now, hardcoded
	BLOCKII.RBRFG.x = 171.835;
	BLOCKII.RBRFG.y = 0.0;
	BLOCKII.RBRFG.z = -10678.596;

	BLOCKII.VBRFG.x = -105.876;
	BLOCKII.VBRFG.y = 0.0;
	BLOCKII.VBRFG.z = -1.04;

	BLOCKII.ABRFG.x = 0.6241;
	BLOCKII.ABRFG.y = 0.0;
	BLOCKII.ABRFG.z = -9.1044;

	BLOCKII.VBRFG_star = -18.72;
	BLOCKII.ABRFG_star = -54.6264;
	BLOCKII.JBRFG_star = -0.01882677;

	BLOCKII.TCGFBRAK = 30.0;
	BLOCKII.TCGIBRAK = 900.0;

	BLOCKII.RAPFG.x = 111.085;
	BLOCKII.RAPFG.y = 0.0;
	BLOCKII.RAPFG.z = -26.794;

	BLOCKII.VAPFG.x = -4.993;
	BLOCKII.VAPFG.y = 0.0;
	BLOCKII.VAPFG.z = 0.248;

	BLOCKII.AAPFG.x = -0.2624;
	BLOCKII.AAPFG.y = 0.0;
	BLOCKII.AAPFG.z = -0.512;

	BLOCKII.VAPFG_star = 4.464;
	BLOCKII.AAPFG_star = -3.072;
	BLOCKII.JAPFG_star = 0.00180772;

	BLOCKII.TCGFAPPR = 30.0;
	BLOCKII.TCGIAPPR = 200.0;

	DescentConstants11_13();

	//TAUVERT
	DoubleToBuffer(10.0*100.0, 14, iTemp, iTemp2);
	SaveEMEM(02516, iTemp);
	SaveEMEM(02517, iTemp2);

	//ABTCOF
	DoubleToBuffer(0.11088229e-6*FT2M/pow(100, 4), -44, iTemp, iTemp2);
	SaveEMEM(02550, iTemp);
	SaveEMEM(02551, iTemp2);

	DoubleToBuffer(-0.58743323e-3*FT2M/pow(100, 3), -27, iTemp, iTemp2);
	SaveEMEM(02552, iTemp);
	SaveEMEM(02553, iTemp2);

	DoubleToBuffer(0.45897234e-1*FT2M/pow(100, 2), -10, iTemp, iTemp2);
	SaveEMEM(02554, iTemp);
	SaveEMEM(02555, iTemp2);

	DoubleToBuffer(5650.9755*FT2M/100.0, 7, iTemp, iTemp2);
	SaveEMEM(02556, iTemp);
	SaveEMEM(02557, iTemp2);

	DoubleToBuffer(0.58577412e-6*FT2M/pow(100, 4), -44, iTemp, iTemp2);
	SaveEMEM(02560, iTemp);
	SaveEMEM(02561, iTemp2);

	DoubleToBuffer(-0.90488751e-3*FT2M/pow(100, 3), -27, iTemp, iTemp2);
	SaveEMEM(02562, iTemp);
	SaveEMEM(02563, iTemp2);

	DoubleToBuffer(0.93189392e-1*FT2M/pow(100, 2), -10, iTemp, iTemp2);
	SaveEMEM(02564, iTemp);
	SaveEMEM(02565, iTemp2);

	DoubleToBuffer(5648.901*FT2M/100.0, 7, iTemp, iTemp2);
	SaveEMEM(02566, iTemp);
	SaveEMEM(02567, iTemp2);

	//VMIN
	DoubleToBuffer(5514.41*FT2M / 100.0, 7, iTemp, iTemp2);
	SaveEMEM(02570, iTemp);
	SaveEMEM(02571, iTemp2);

	//YLIM
	DoubleToBuffer(8.2*1852.0, 24, iTemp, iTemp2);
	SaveEMEM(02572, iTemp);
	SaveEMEM(02573, iTemp2);

	//ABTRDOT
	DoubleToBuffer(19.5*0.3048 / 100.0, 7, iTemp, iTemp2);
	SaveEMEM(02574, iTemp);
	SaveEMEM(02575, iTemp2);

	//COSTHET1
	SaveEMEM(02576, 0);
	SaveEMEM(02577, 0);

	//COSTHET2
	SaveEMEM(02600, 06733);
	SaveEMEM(02601, 07535);

	//TNEWA
	SaveEMEM(03431, 0);
	SaveEMEM(03432, 017500);
}

void AGCPadloadGenerator::Luminary116Padload()
{
	LGCDefaults(false);
	Luminary099_116_Defaults();

	//FLAGWRD3
	SaveEMEM(077, 012000);

	//MASS
	DoubleToBuffer(BLOCKII.TotalMass*LBS2KG, 16, iTemp, iTemp2);
	SaveEMEM(01243, iTemp);
	SaveEMEM(01244, 0);

	//LEMMASS
	iTemp = SingleToBuffer(BLOCKII.LMMass*LBS2KG, 16);
	SaveEMEM(01326, iTemp);
	//CSMMASS
	iTemp = SingleToBuffer(BLOCKII.CSMMass*LBS2KG, 16);
	SaveEMEM(01327, iTemp);

	R2Model(01347);

	//RADSKAL
	SaveEMEM(01351, 0);
	SaveEMEM(01352, 0);
	//SKALSKAL
	SaveEMEM(01353, 0);

	//REFSMMAT
	SaveEMEM(01733, 061727);
	SaveEMEM(01734, 040512);
	SaveEMEM(01735, 075331);
	SaveEMEM(01736, 057404);
	SaveEMEM(01737, 07114);
	SaveEMEM(01740, 027444);
	SaveEMEM(01741, 07225);
	SaveEMEM(01742, 032063);
	SaveEMEM(01743, 067137);
	SaveEMEM(01744, 063015);
	SaveEMEM(01745, 013137);
	SaveEMEM(01746, 023771);
	SaveEMEM(01747, 02106);
	SaveEMEM(01750, 022336);
	SaveEMEM(01751, 015064);
	SaveEMEM(01752, 032774);
	SaveEMEM(01753, 010733);
	SaveEMEM(01754, 034410);

	//For now, hardcoded
	BLOCKII.RBRFG.x = -3562.05;
	BLOCKII.RBRFG.y = 0.0;
	BLOCKII.RBRFG.z = -13705.71;

	BLOCKII.VBRFG.x = -186.90305;
	BLOCKII.VBRFG.y = 0.0;
	BLOCKII.VBRFG.z = -98.73819;

	BLOCKII.ABRFG.x = -0.4502495;
	BLOCKII.ABRFG.y = 0.0;
	BLOCKII.ABRFG.z = -9.5150975;

	BLOCKII.VBRFG_star = -1777.28742;
	BLOCKII.ABRFG_star = -57.090585;
	BLOCKII.JBRFG_star = -1.47427e-2;

	BLOCKII.TCGFBRAK = 30.0;
	BLOCKII.TCGIBRAK = 900.0;

	BLOCKII.RAPFG.x = 82.9275;
	BLOCKII.RAPFG.y = 0.0;
	BLOCKII.RAPFG.z = -20.1605;

	BLOCKII.VAPFG.x = -0.319;
	BLOCKII.VAPFG.y = 0.0;
	BLOCKII.VAPFG.z = 0.31233;

	BLOCKII.AAPFG.x = 0.29982;
	BLOCKII.AAPFG.y = 0.0;
	BLOCKII.AAPFG.z = -0.40165;

	BLOCKII.VAPFG_star = 5.62194;
	BLOCKII.AAPFG_star = -2.4099;
	BLOCKII.JAPFG_star = 0.0376954;

	BLOCKII.TCGFAPPR = 6.0;
	BLOCKII.TCGIAPPR = 200.0;

	DescentConstants11_13();

	//TAUVERT
	DoubleToBuffer(8.0*100.0, 14, iTemp, iTemp2);
	SaveEMEM(02516, iTemp);
	SaveEMEM(02517, iTemp2);

	//J1PARM
	DoubleToBuffer(BLOCKII.J1PARM*FT2M, 23, iTemp, iTemp2);
	SaveEMEM(02550, iTemp);
	SaveEMEM(02551, iTemp2);
	//K1PARM
	DoubleToBuffer(BLOCKII.K1PARM*FT2M*PI2, 23, iTemp, iTemp2);
	SaveEMEM(02552, iTemp);
	SaveEMEM(02553, iTemp2);
	//J2PARM
	DoubleToBuffer(BLOCKII.J2PARM*FT2M, 23, iTemp, iTemp2);
	SaveEMEM(02554, iTemp);
	SaveEMEM(02555, iTemp2);
	//K2PARM
	DoubleToBuffer(BLOCKII.K2PARM*FT2M*PI2, 23, iTemp, iTemp2);
	SaveEMEM(02556, iTemp);
	SaveEMEM(02557, iTemp2);
	//THETCRIT
	DoubleToBuffer(BLOCKII.THETCRIT / 360.0, 0, iTemp, iTemp2);
	SaveEMEM(02560, iTemp);
	SaveEMEM(02561, iTemp2);
	//RAMIN
	DoubleToBuffer(BLOCKII.RAMIN*0.3048, 24, iTemp, iTemp2);
	SaveEMEM(02562, iTemp);
	SaveEMEM(02563, iTemp2);
	//YLIM
	DoubleToBuffer(8.2*1852.0, 24, iTemp, iTemp2);
	SaveEMEM(02564, iTemp);
	SaveEMEM(02565, iTemp2);
	//ABTRDOT
	DoubleToBuffer(19.5*0.3048 / 100.0, 7, iTemp, iTemp2);
	SaveEMEM(02566, iTemp);
	SaveEMEM(02567, iTemp2);
	//COSTHET1
	SaveEMEM(02570, 0);
	SaveEMEM(02571, 0);
	//COSTHET2
	SaveEMEM(02572, 06733);
	SaveEMEM(02573, 07535);

	//DLAND
	SaveEMEM(02634, 0);
	SaveEMEM(02635, 0);
	SaveEMEM(02636, 0);
	SaveEMEM(02637, 0);
	SaveEMEM(02640, 0);
	SaveEMEM(02641, 0);

	//TNEWA
	SaveEMEM(03431, 020000);
	SaveEMEM(03432, 0);
}

void AGCPadloadGenerator::Luminary131Padload()
{
	LGCDefaults(false);

	//FLAGWRD3
	SaveEMEM(077, 012000);
	//FLAGWRD8
	SaveEMEM(0104, 06000);
	//FLAGWRD10
	SaveEMEM(0106, 0);

	//MASS
	DoubleToBuffer(BLOCKII.TotalMass*LBS2KG, 16, iTemp, iTemp2);
	SaveEMEM(01243, iTemp);
	SaveEMEM(01244, 0);

	//LEMMASS
	iTemp = SingleToBuffer(BLOCKII.LMMass*LBS2KG, 16);
	SaveEMEM(01326, iTemp);
	//CSMMASS
	iTemp = SingleToBuffer(BLOCKII.CSMMass*LBS2KG, 16);
	SaveEMEM(01327, iTemp);

	R2Model(01347);

	//RADSKAL
	SaveEMEM(01351, 0);
	SaveEMEM(01352, 0);
	//SKALSKAL
	SaveEMEM(01353, 0);

	//ELBIAS
	SaveEMEM(01356, 0);

	//TETCSM
	SaveEMEM(01570, 037777);
	//TETLEM
	SaveEMEM(01642, 037777);

	//REFSMMAT
	SaveEMEM(01733, 011227);
	SaveEMEM(01734, 03621);
	SaveEMEM(01735, 064303);
	SaveEMEM(01736, 076432);
	SaveEMEM(01737, 072072);
	SaveEMEM(01740, 064034);
	SaveEMEM(01741, 077664);
	SaveEMEM(01742, 050201);
	SaveEMEM(01743, 070531);
	SaveEMEM(01744, 072041);
	SaveEMEM(01745, 016162);
	SaveEMEM(01746, 031571);
	SaveEMEM(01747, 062764);
	SaveEMEM(01750, 063727);
	SaveEMEM(01751, 067624);
	SaveEMEM(01752, 076002);
	SaveEMEM(01753, 073506);
	SaveEMEM(01754, 053045);

	//RANGEVAR
	DoubleToBuffer(1.111111111e-5, -12, iTemp, iTemp2);
	SaveEMEM(01770, iTemp);
	SaveEMEM(01771, iTemp2);
	//RATEVAR
	DoubleToBuffer(1.877777000e-5, -12, iTemp, iTemp2);
	SaveEMEM(01772, iTemp);
	SaveEMEM(01773, iTemp2);
	//RVARMIN
	iTemp = SingleToBuffer(66.0, 12);
	SaveEMEM(01774, iTemp);
	//VVARMIN
	iTemp = SingleToBuffer(1.7445e-6, -12);
	SaveEMEM(01775, iTemp);

	//WSURFPOS
	SaveEMEM(02006, 0);
	//WSURFVEL
	SaveEMEM(02007, 0);

	//SHAFTVAR
	iTemp = SingleToBuffer(BLOCKII.SHAFTVAR*1e-6, -12);
	SaveEMEM(02010, iTemp);
	//TRUNVAR
	iTemp = SingleToBuffer(BLOCKII.TRUNVAR*1e-6, -12);
	SaveEMEM(02011, iTemp);

	//AGSK
	DoubleToBuffer(BLOCKII.AGSK*3600.0*100.0, 28, iTemp, iTemp2);
	SaveEMEM(02020, iTemp);
	SaveEMEM(02021, iTemp2);

	//RLS
	RLS = r_from_latlong(LSLat*RAD, LSLng*RAD, LSAlt, BODY_MOON, 1);
	DoubleToBuffer(RLS.x, 27, iTemp, iTemp2);
	SaveEMEM(02022, iTemp);
	SaveEMEM(02023, iTemp2);
	DoubleToBuffer(RLS.y, 27, iTemp, iTemp2);
	SaveEMEM(02024, iTemp);
	SaveEMEM(02025, iTemp2);
	DoubleToBuffer(RLS.z, 27, iTemp, iTemp2);
	SaveEMEM(02026, iTemp);
	SaveEMEM(02027, iTemp2);

	//For now, hardcoded
	BLOCKII.RBRFG.x = -3562.05;
	BLOCKII.RBRFG.y = 0.0;
	BLOCKII.RBRFG.z = -13705.71;

	BLOCKII.VBRFG.x = -186.90305;
	BLOCKII.VBRFG.y = 0.0;
	BLOCKII.VBRFG.z = -98.73819;

	BLOCKII.ABRFG.x = -0.4502495;
	BLOCKII.ABRFG.y = 0.0;
	BLOCKII.ABRFG.z = -9.5150975;

	BLOCKII.VBRFG_star = -1777.28742;
	BLOCKII.ABRFG_star = -57.090585;
	BLOCKII.JBRFG_star = -1.4742736e-2;

	BLOCKII.TCGFBRAK = 30.0;
	BLOCKII.TCGIBRAK = 900.0;

	BLOCKII.RAPFG.x = 82.9275;
	BLOCKII.RAPFG.y = 0.0;
	BLOCKII.RAPFG.z = -20.1605;

	BLOCKII.RAPFG.x = 82.9275;
	BLOCKII.RAPFG.y = 0.0;
	BLOCKII.RAPFG.z = -20.1605;

	BLOCKII.VAPFG.x = -0.319;
	BLOCKII.VAPFG.y = 0.0;
	BLOCKII.VAPFG.z = 0.31233;

	BLOCKII.AAPFG.x = 0.29982;
	BLOCKII.AAPFG.y = 0.0;
	BLOCKII.AAPFG.z = -0.40165;

	BLOCKII.VAPFG_star = 5.62194;
	BLOCKII.AAPFG_star = -2.4099;
	BLOCKII.JAPFG_star = 0.03769542;

	BLOCKII.TCGFAPPR = 6.0;
	BLOCKII.TCGIAPPR = 200.0;

	DescentConstants11_13();

	//TAUHZ
	iTemp = SingleToBuffer(5.0*100.0, 11);
	SaveEMEM(02510, iTemp);
	//QHZ
	SaveEMEM(02511, 014632);
	//AHZLIM
	SaveEMEM(02512, 017);
	//TOOFEW
	SaveEMEM(02513, 03);
	//HLROFF
	DoubleToBuffer(15.25, 24, iTemp, iTemp2);
	SaveEMEM(02514, iTemp);
	SaveEMEM(02515, iTemp2);
	//2LATE466
	DoubleToBuffer(1.50*100.0, 28, iTemp, iTemp2);
	SaveEMEM(02516, iTemp);
	SaveEMEM(02517, iTemp2);

	//J1PARM
	DoubleToBuffer(BLOCKII.J1PARM*FT2M, 23, iTemp, iTemp2);
	SaveEMEM(02550, iTemp);
	SaveEMEM(02551, iTemp2);
	//K1PARM
	DoubleToBuffer(BLOCKII.K1PARM *FT2M*PI2, 23, iTemp, iTemp2);
	SaveEMEM(02552, iTemp);
	SaveEMEM(02553, iTemp2);
	//J2PARM
	DoubleToBuffer(BLOCKII.J2PARM*FT2M, 23, iTemp, iTemp2);
	SaveEMEM(02554, iTemp);
	SaveEMEM(02555, iTemp2);
	//K2PARM
	DoubleToBuffer(BLOCKII.K2PARM *FT2M*PI2, 23, iTemp, iTemp2);
	SaveEMEM(02556, iTemp);
	SaveEMEM(02557, iTemp2);
	//THETCRIT
	DoubleToBuffer(BLOCKII.THETCRIT / 360.0, 0, iTemp, iTemp2);
	SaveEMEM(02560, iTemp);
	SaveEMEM(02561, iTemp2);
	//RAMIN
	DoubleToBuffer(BLOCKII.RAMIN*0.3048, 24, iTemp, iTemp2);
	SaveEMEM(02562, iTemp);
	SaveEMEM(02563, iTemp2);
	//YLIM
	DoubleToBuffer(8.2*1852.0, 24, iTemp, iTemp2);
	SaveEMEM(02564, iTemp);
	SaveEMEM(02565, iTemp2);
	//ABTRDOT
	DoubleToBuffer(19.5*0.3048/100.0, 7, iTemp, iTemp2);
	SaveEMEM(02566, iTemp);
	SaveEMEM(02567, iTemp2);
	//COSTHET1
	SaveEMEM(02570, 0);
	SaveEMEM(02571, 0);
	//COSTHET2
	SaveEMEM(02572, 06733);
	SaveEMEM(02573, 07535);

	//DLAND
	SaveEMEM(02634, 0);
	SaveEMEM(02635, 0);
	SaveEMEM(02636, 0);
	SaveEMEM(02637, 0);
	SaveEMEM(02640, 0);
	SaveEMEM(02641, 0);

	//DOWNTORK
	SaveEMEM(03113, 0);
	SaveEMEM(03114, 0);
	SaveEMEM(03115, 0);
	SaveEMEM(03116, 0);
	SaveEMEM(03117, 0);
	SaveEMEM(03120, 0);

	//VELBIAS
	DoubleToBuffer(2.5*0.3048/100.0, 6, iTemp, iTemp2);
	SaveEMEM(03371, iTemp);
	SaveEMEM(03372, iTemp2);

	//AZBIAS
	iTemp = SingleToBuffer(0.0 / 360.0, -1);
	SaveEMEM(03373, iTemp);

	//AOTAZ
	iTemp = SingleToBuffer(-60.0 / 360.0, -1, true);
	SaveEMEM(03404, iTemp);
	iTemp = SingleToBuffer(0.0 / 360.0, -1, true);
	SaveEMEM(03405, iTemp);
	iTemp = SingleToBuffer(60.0 / 360.0, -1, true);
	SaveEMEM(03406, iTemp);
	iTemp = SingleToBuffer(120.0 / 360.0, -1, true);
	SaveEMEM(03407, iTemp);
	iTemp = SingleToBuffer(-180.0 / 360.0, -1, true);
	SaveEMEM(03410, iTemp);
	iTemp = SingleToBuffer(-120.0 / 360.0, -1, true);
	SaveEMEM(03411, iTemp);
	//AOTEL
	iTemp = SingleToBuffer(45.0 / 360.0, -1);
	SaveEMEM(03412, iTemp);
	iTemp = SingleToBuffer(45.0 / 360.0, -1);
	SaveEMEM(03413, iTemp);
	iTemp = SingleToBuffer(45.0 / 360.0, -1);
	SaveEMEM(03414, iTemp);
	iTemp = SingleToBuffer(45.0 / 360.0, -1);
	SaveEMEM(03415, iTemp);
	iTemp = SingleToBuffer(45.0 / 360.0, -1);
	SaveEMEM(03416, iTemp);
	iTemp = SingleToBuffer(45.0 / 360.0, -1);
	SaveEMEM(03417, iTemp);

	//LRHMAX
	iTemp = SingleToBuffer(50000.0*0.3048, 14);
	SaveEMEM(03420, iTemp);
	//LRWH
	iTemp = SingleToBuffer(0.35, 0);
	SaveEMEM(03421, iTemp);
	//ZOOMTIME
	iTemp = SingleToBuffer(26.0*100.0, 14);
	SaveEMEM(03422, iTemp);
	//TENDBRAK
	iTemp = SingleToBuffer(62.0*100.0, 17);
	SaveEMEM(03423, iTemp);
	//TENDAPPR
	iTemp = SingleToBuffer(12.0*100.0, 17);
	SaveEMEM(03424, iTemp);

	//LEADTIME
	iTemp = SingleToBuffer(-2.2*100.0, 17);
	SaveEMEM(03426, iTemp);
	//LEADTIME
	iTemp = SingleToBuffer(62.0*100.0, 17);
	SaveEMEM(03427, iTemp);
	//RPCRTQSW
	iTemp = SingleToBuffer(-1.0, 1);
	SaveEMEM(03430, iTemp);
	//TNEWA
	SaveEMEM(03431, 020000);
	SaveEMEM(03432, 0);
}

void AGCPadloadGenerator::Luminary099_116_Defaults()
{
	//Same addresses and values for Luminary 99 and 116

	//FLAGWRD8
	SaveEMEM(0104, 0);
	//FLAGWRD10
	SaveEMEM(0106, 0);

	//TETCSM
	SaveEMEM(01570, 037777);
	//TETLEM
	SaveEMEM(01642, 037777);

	//RANGEVAR
	DoubleToBuffer(1.111111111e-5, -12, iTemp, iTemp2);
	SaveEMEM(01770, iTemp);
	SaveEMEM(01771, iTemp2);
	//RATEVAR
	DoubleToBuffer(1.877777000e-5, -12, iTemp, iTemp2);
	SaveEMEM(01772, iTemp);
	SaveEMEM(01773, iTemp2);
	//RVARMIN
	iTemp = SingleToBuffer(66.0, 12);
	SaveEMEM(01774, iTemp);
	//VVARMIN
	iTemp = SingleToBuffer(1.7445e-6, -12);
	SaveEMEM(01775, iTemp);

	//WSURFPOS
	iTemp = SingleToBuffer(BLOCKII.WSURFPOS*FT2M, 14);
	SaveEMEM(02006, iTemp);
	//WSURFVEL
	iTemp = SingleToBuffer(BLOCKII.WSURFVEL*FT2M / 100.0, 0);
	SaveEMEM(02007, iTemp);

	//SHAFTVAR
	iTemp = SingleToBuffer(BLOCKII.SHAFTVAR*1e-6, -12);
	SaveEMEM(02010, iTemp);
	//TRUNVAR
	iTemp = SingleToBuffer(BLOCKII.TRUNVAR*1e-6, -12);
	SaveEMEM(02011, iTemp);

	//AGSK
	DoubleToBuffer(BLOCKII.AGSK*3600.0*100.0, 28, iTemp, iTemp2);
	SaveEMEM(02020, iTemp);
	SaveEMEM(02021, iTemp2);

	//RLS
	RLS = r_from_latlong(LSLat*RAD, LSLng*RAD, LSAlt, BODY_MOON, 1);
	DoubleToBuffer(RLS.x, 27, iTemp, iTemp2);
	SaveEMEM(02022, iTemp);
	SaveEMEM(02023, iTemp2);
	DoubleToBuffer(RLS.y, 27, iTemp, iTemp2);
	SaveEMEM(02024, iTemp);
	SaveEMEM(02025, iTemp2);
	DoubleToBuffer(RLS.z, 27, iTemp, iTemp2);
	SaveEMEM(02026, iTemp);
	SaveEMEM(02027, iTemp2);

	//V2FG
	DoubleToBuffer(-3.0*FT2M / 100.0, 10, iTemp, iTemp2);
	SaveEMEM(02510, iTemp);
	SaveEMEM(02511, iTemp2);
	DoubleToBuffer(0.0*FT2M / 100.0, 10, iTemp, iTemp2);
	SaveEMEM(02512, iTemp);
	SaveEMEM(02513, iTemp2);
	DoubleToBuffer(0.0*FT2M / 100.0, 10, iTemp, iTemp2);
	SaveEMEM(02514, iTemp);
	SaveEMEM(02515, iTemp2);

	//DOWNTORK
	SaveEMEM(03113, 0);
	SaveEMEM(03114, 0);
	SaveEMEM(03115, 0);
	SaveEMEM(03116, 0);
	SaveEMEM(03117, 0);
	SaveEMEM(03120, 0);

	//AOTAZ
	iTemp = SingleToBuffer(-60.0 / 360.0, -1, true);
	SaveEMEM(03404, iTemp);
	iTemp = SingleToBuffer(0.0 / 360.0, -1, true);
	SaveEMEM(03405, iTemp);
	iTemp = SingleToBuffer(60.0 / 360.0, -1, true);
	SaveEMEM(03406, iTemp);
	iTemp = SingleToBuffer(120.0 / 360.0, -1, true);
	SaveEMEM(03407, iTemp);
	iTemp = SingleToBuffer(-180.0 / 360.0, -1, true);
	SaveEMEM(03410, iTemp);
	iTemp = SingleToBuffer(-120.0 / 360.0, -1, true);
	SaveEMEM(03411, iTemp);
	//AOTEL
	iTemp = SingleToBuffer(45.0 / 360.0, -1);
	SaveEMEM(03412, iTemp);
	iTemp = SingleToBuffer(45.0 / 360.0, -1);
	SaveEMEM(03413, iTemp);
	iTemp = SingleToBuffer(45.0 / 360.0, -1);
	SaveEMEM(03414, iTemp);
	iTemp = SingleToBuffer(45.0 / 360.0, -1);
	SaveEMEM(03415, iTemp);
	iTemp = SingleToBuffer(45.0 / 360.0, -1);
	SaveEMEM(03416, iTemp);
	iTemp = SingleToBuffer(45.0 / 360.0, -1);
	SaveEMEM(03417, iTemp);

	//LRHMAX
	iTemp = SingleToBuffer(50000.0*0.3048, 14);
	SaveEMEM(03420, iTemp);
	//LRWH
	iTemp = SingleToBuffer(0.35, 0);
	SaveEMEM(03421, iTemp);
	//ZOOMTIME
	iTemp = SingleToBuffer(26.0*100.0, 14);
	SaveEMEM(03422, iTemp);
	//TENDBRAK
	iTemp = SingleToBuffer(62.0*100.0, 17);
	SaveEMEM(03423, iTemp);
	//TENDAPPR
	iTemp = SingleToBuffer(12.0*100.0, 17);
	SaveEMEM(03424, iTemp);

	//LEADTIME
	iTemp = SingleToBuffer(-2.2*100.0, 17);
	SaveEMEM(03426, iTemp);
	//LEADTIME
	iTemp = SingleToBuffer(62.0*100.0, 17);
	SaveEMEM(03427, iTemp);
	//RPCRTQSW
	iTemp = SingleToBuffer(-1.0, 1);
	SaveEMEM(03430, iTemp);
}

void AGCPadloadGenerator::DescentConstants11_13()
{
	//Same addresses from Apollo 11 to 13

	//TLAND
	DoubleToBuffer(BLOCKII.TLAND*3600.0*100.0, 28, iTemp, iTemp2);
	SaveEMEM(02400, iTemp);
	SaveEMEM(02401, iTemp2);

	//RBRFG
	DoubleToBuffer(BLOCKII.RBRFG.x*FT2M, 24, iTemp, iTemp2);
	SaveEMEM(02402, iTemp);
	SaveEMEM(02403, iTemp2);
	DoubleToBuffer(BLOCKII.RBRFG.y*FT2M, 24, iTemp, iTemp2);
	SaveEMEM(02404, iTemp);
	SaveEMEM(02405, iTemp2);
	DoubleToBuffer(BLOCKII.RBRFG.z*FT2M, 24, iTemp, iTemp2);
	SaveEMEM(02406, iTemp);
	SaveEMEM(02407, iTemp2);

	//VBRFG
	DoubleToBuffer(BLOCKII.VBRFG.x*FT2M / 100.0, 10, iTemp, iTemp2);
	SaveEMEM(02410, iTemp);
	SaveEMEM(02411, iTemp2);
	DoubleToBuffer(BLOCKII.VBRFG.y*FT2M / 100.0, 10, iTemp, iTemp2);
	SaveEMEM(02412, iTemp);
	SaveEMEM(02413, iTemp2);
	DoubleToBuffer(BLOCKII.VBRFG.z*FT2M / 100.0, 10, iTemp, iTemp2);
	SaveEMEM(02414, iTemp);
	SaveEMEM(02415, iTemp2);

	//ABRFG
	DoubleToBuffer(BLOCKII.ABRFG.x*FT2M / 10000.0, -4, iTemp, iTemp2);
	SaveEMEM(02416, iTemp);
	SaveEMEM(02417, iTemp2);
	DoubleToBuffer(BLOCKII.ABRFG.y*FT2M / 10000.0, -4, iTemp, iTemp2);
	SaveEMEM(02420, iTemp);
	SaveEMEM(02421, iTemp2);
	DoubleToBuffer(BLOCKII.ABRFG.z*FT2M / 10000.0, -4, iTemp, iTemp2);
	SaveEMEM(02422, iTemp);
	SaveEMEM(02423, iTemp2);

	//VBRFG*
	DoubleToBuffer(BLOCKII.VBRFG_star*FT2M / 100.0, 13, iTemp, iTemp2);
	SaveEMEM(02424, iTemp);
	SaveEMEM(02425, iTemp2);
	//ABRFG*
	DoubleToBuffer(BLOCKII.ABRFG_star*FT2M / 10000.0, -4, iTemp, iTemp2);
	SaveEMEM(02426, iTemp);
	SaveEMEM(02427, iTemp2);
	//JBRFG*
	DoubleToBuffer(BLOCKII.JBRFG_star*FT2M / 1000000.0, -21, iTemp, iTemp2);
	SaveEMEM(02430, iTemp);
	SaveEMEM(02431, iTemp2);

	//GAINBRAK
	SaveEMEM(02432, 037777);
	SaveEMEM(02433, 037777);
	//TCGFBRAK
	iTemp = SingleToBuffer(BLOCKII.TCGFBRAK*100.0, 17);
	SaveEMEM(02434, iTemp);
	//TCGIBRAK
	iTemp = SingleToBuffer(BLOCKII.TCGIBRAK*100.0, 17);
	SaveEMEM(02435, iTemp);

	//RAPFG
	DoubleToBuffer(BLOCKII.RAPFG.x*FT2M, 24, iTemp, iTemp2);
	SaveEMEM(02436, iTemp);
	SaveEMEM(02437, iTemp2);
	DoubleToBuffer(BLOCKII.RAPFG.y*FT2M, 24, iTemp, iTemp2);
	SaveEMEM(02440, iTemp);
	SaveEMEM(02441, iTemp2);
	DoubleToBuffer(BLOCKII.RAPFG.z*FT2M, 24, iTemp, iTemp2);
	SaveEMEM(02442, iTemp);
	SaveEMEM(02443, iTemp2);

	//VARFG
	DoubleToBuffer(BLOCKII.VAPFG.x*FT2M / 100.0, 10, iTemp, iTemp2);
	SaveEMEM(02444, iTemp);
	SaveEMEM(02445, iTemp2);
	DoubleToBuffer(BLOCKII.VAPFG.y*FT2M / 100.0, 10, iTemp, iTemp2);
	SaveEMEM(02446, iTemp);
	SaveEMEM(02447, iTemp2);
	DoubleToBuffer(BLOCKII.VAPFG.z*FT2M / 100.0, 10, iTemp, iTemp2);
	SaveEMEM(02450, iTemp);
	SaveEMEM(02451, iTemp2);

	//AARFG
	DoubleToBuffer(BLOCKII.AAPFG.x*FT2M / 10000.0, -4, iTemp, iTemp2);
	SaveEMEM(02452, iTemp);
	SaveEMEM(02453, iTemp2);
	DoubleToBuffer(BLOCKII.AAPFG.y*FT2M / 10000.0, -4, iTemp, iTemp2);
	SaveEMEM(02454, iTemp);
	SaveEMEM(02455, iTemp2);
	DoubleToBuffer(BLOCKII.AAPFG.z*FT2M / 10000.0, -4, iTemp, iTemp2);
	SaveEMEM(02456, iTemp);
	SaveEMEM(02457, iTemp2);

	//VARFG*
	DoubleToBuffer(BLOCKII.VAPFG_star*FT2M / 100.0, 13, iTemp, iTemp2);
	SaveEMEM(02460, iTemp);
	SaveEMEM(02461, iTemp2);
	//AARFG*
	DoubleToBuffer(BLOCKII.AAPFG_star*FT2M / 10000.0, -4, iTemp, iTemp2);
	SaveEMEM(02462, iTemp);
	SaveEMEM(02463, iTemp2);
	//JARFG*
	DoubleToBuffer(BLOCKII.JAPFG_star*FT2M / 1000000.0, -21, iTemp, iTemp2);
	SaveEMEM(02464, iTemp);
	SaveEMEM(02465, iTemp2);

	//GAINAPPR
	SaveEMEM(02466, 0);
	SaveEMEM(02467, 0);
	//TCGFAPPR
	iTemp = SingleToBuffer(BLOCKII.TCGFAPPR*100.0, 17);
	SaveEMEM(02470, iTemp);
	//TCGIAPPR
	iTemp = SingleToBuffer(BLOCKII.TCGIAPPR*100.0, 17);
	SaveEMEM(02471, iTemp);

	//VIGN
	DoubleToBuffer(BLOCKII.VIGN*FT2M / 100.0, 10, iTemp, iTemp2);
	SaveEMEM(02472, iTemp);
	SaveEMEM(02473, iTemp2);
	//RIGNX
	DoubleToBuffer(BLOCKII.RIGNX*0.3048, 24, iTemp, iTemp2);
	SaveEMEM(02474, iTemp);
	SaveEMEM(02475, iTemp2);
	//RIGNZ
	DoubleToBuffer(BLOCKII.RIGNZ*0.3048, 24, iTemp, iTemp2);
	SaveEMEM(02476, iTemp);
	SaveEMEM(02477, iTemp2);
	//KIGNX/B4
	DoubleToBuffer(BLOCKII.KIGNXB4, 4, iTemp, iTemp2);
	SaveEMEM(02500, iTemp);
	SaveEMEM(02501, iTemp2);
	//KIGNY/B8
	DoubleToBuffer(BLOCKII.KIGNYB8 / 0.3048, -16, iTemp, iTemp2);
	SaveEMEM(02502, iTemp);
	SaveEMEM(02503, iTemp2);
	//KIGNV/B4
	DoubleToBuffer(BLOCKII.KIGNVB4*100.0, 18, iTemp, iTemp2);
	SaveEMEM(02504, iTemp);
	SaveEMEM(02505, iTemp2);
	//LOWCRIT
	SaveEMEM(02506, 04114);
	//HIGHCRIT
	SaveEMEM(02507, 04454);

	//DELQFIX
	DoubleToBuffer(200.0*0.3048, 24, iTemp, iTemp2);
	SaveEMEM(02520, iTemp);
	SaveEMEM(02521, iTemp2);
	//LRALPHA
	iTemp = SingleToBuffer(6.0 / 360.0, -1);
	SaveEMEM(02522, iTemp);
	//LRBETA1
	iTemp = SingleToBuffer(24.0 / 360.0, -1);
	SaveEMEM(02523, iTemp);
	//LRALPHA2
	iTemp = SingleToBuffer(6.0 / 360.0, -1);
	SaveEMEM(02524, iTemp);
	//LRBETA2
	iTemp = SingleToBuffer(0.0 / 360.0, -1);
	SaveEMEM(02525, iTemp);
	//LRVMAX
	iTemp = SingleToBuffer(2000.0*0.3048 / 100.0, 7);
	SaveEMEM(02526, iTemp);
	//LRVF
	iTemp = SingleToBuffer(200.0*0.3048 / 100.0, 7);
	SaveEMEM(02527, iTemp);
	//LRWVZ
	iTemp = SingleToBuffer(0.3, 0);
	SaveEMEM(02530, iTemp);
	//LRWVY
	iTemp = SingleToBuffer(0.3, 0);
	SaveEMEM(02531, iTemp);
	//LRWVX
	iTemp = SingleToBuffer(0.3, 0);
	SaveEMEM(02532, iTemp);
	//LRWVFZ
	iTemp = SingleToBuffer(0.2, 0);
	SaveEMEM(02533, iTemp);
	//LRWVFY
	iTemp = SingleToBuffer(0.2, 0);
	SaveEMEM(02534, iTemp);
	//LRWVFX
	iTemp = SingleToBuffer(0.2, 0);
	SaveEMEM(02535, iTemp);
	//LRWVFF
	iTemp = SingleToBuffer(0.1, 0);
	SaveEMEM(02536, iTemp);
	//RODSCALE
	iTemp = SingleToBuffer(1.0*0.3048 / 100.0, -7);
	SaveEMEM(02537, iTemp);
	//TAUROD
	DoubleToBuffer(1.5*100.0, 9, iTemp, iTemp2);
	SaveEMEM(02540, iTemp);
	SaveEMEM(02541, iTemp2);
	//LAG/TAU
	DoubleToBuffer(0.413333, 0, iTemp, iTemp2);
	SaveEMEM(02542, iTemp);
	SaveEMEM(02543, iTemp2);
	//MINFORCE
	SaveEMEM(02544, 01);
	SaveEMEM(02545, 027631);
	//MAXFORCE
	SaveEMEM(02546, 013);
	SaveEMEM(02547, 06551);

	//IGNAOSQ
	iTemp = SingleToBuffer(BLOCKII.IGNAOSQ / 360.0, -2);
	SaveEMEM(03012, iTemp);
	//IGNAOSR
	iTemp = SingleToBuffer(BLOCKII.IGNAOSR / 360.0, -2);
	SaveEMEM(03013, iTemp);

	//DELTTFAP
	iTemp = SingleToBuffer(BLOCKII.DELTTFAP*100.0, 17);
	SaveEMEM(03425, iTemp);
}

void AGCPadloadGenerator::Luminary178Padload()
{
	LGCDefaults(false);

	//REFSMMAT
	SaveEMEM(01731, 067061);
	SaveEMEM(01732, 060375);
	SaveEMEM(01733, 064033);
	SaveEMEM(01734, 051650);
	SaveEMEM(01735, 072126);
	SaveEMEM(01736, 045456);
	SaveEMEM(01737, 02656);
	SaveEMEM(01740, 014434);
	SaveEMEM(01741, 067336);
	SaveEMEM(01742, 043072);
	SaveEMEM(01743, 015154);
	SaveEMEM(01744, 015512);
	SaveEMEM(01745, 063006);
	SaveEMEM(01746, 064677);
	SaveEMEM(01747, 06242);
	SaveEMEM(01750, 0131);
	SaveEMEM(01751, 06706);
	SaveEMEM(01752, 020152);

	BLOCKII.RBRFG.x = -1773.725;
	BLOCKII.RAPFG.x = 94.91910105;
	BLOCKII.RBRFG.z = -14488.0269;
	BLOCKII.RAPFG.z = -15.72079987;
	BLOCKII.VBRFG.x = -168.1064567;
	BLOCKII.VAPFG.x = 2.083579987;
	BLOCKII.VBRFG.z = -77.6143668;
	BLOCKII.VAPFG.z = 0.8303187992;
	BLOCKII.ABRFG.x = -0.6472360236;
	BLOCKII.AAPFG.x = 0.5402850066;
	BLOCKII.ABRFG.z = -8.41438189;
	BLOCKII.AAPFG.z = -0.2354229987;
	BLOCKII.VBRFG_star = -1397.058596;
	BLOCKII.VAPFG_star = 14.94573786;
	BLOCKII.ABRFG_star = -50.48629265;
	BLOCKII.AAPFG_star = -1.412537992;
	BLOCKII.JBRFG_star = 0.008257294947;
	BLOCKII.JAPFG_star = 0.04509242126;

	BLOCKII.TCGFBRAK = 30.0;
	BLOCKII.TCGIBRAK = 900.0;
	BLOCKII.TCGFAPPR = 6.0;
	BLOCKII.TCGIAPPR = 200.0;

	BLOCKII.DELQFIX = 500.0;

	DescentConstants14_17();

	//LRWH1
	iTemp = SingleToBuffer(0.35, 0);
	SaveEMEM(03756, iTemp);
}

void AGCPadloadGenerator::Zerlina56Padload()
{
	LGCDefaults(false);

	//FLAGWRD3
	SaveEMEM(077, 012000);
	//FLAGWRD8
	SaveEMEM(0104, 06000);
	//FLAGWRD10
	SaveEMEM(0106, 0);

	//MASS
	DoubleToBuffer(BLOCKII.TotalMass*LBS2KG, 16, iTemp, iTemp2);
	SaveEMEM(01245, iTemp);
	SaveEMEM(01246, 0);

	//LRWH1
	iTemp = SingleToBuffer(0.35, 0);
	SaveEMEM(01315, iTemp);

	//LEMMASS
	iTemp = SingleToBuffer(BLOCKII.LMMass*LBS2KG, 16);
	SaveEMEM(01326, iTemp);
	//CSMMASS
	iTemp = SingleToBuffer(BLOCKII.CSMMass*LBS2KG, 16);
	SaveEMEM(01327, iTemp);

	R2Model(01347);

	//ELBIAS
	SaveEMEM(01353, 0);

	//TOOFEW
	SaveEMEM(01354, 03);

	//TETCSM
	SaveEMEM(01570, 037777);
	SaveEMEM(01571, 037777);
	//TETLEM
	SaveEMEM(01642, 037777);
	SaveEMEM(01643, 037777);

	//REFSMMAT (Just a generic one)
	SaveEMEM(01731, 067061);
	SaveEMEM(01732, 060375);
	SaveEMEM(01733, 064033);
	SaveEMEM(01734, 051650);
	SaveEMEM(01735, 072126);
	SaveEMEM(01736, 045456);
	SaveEMEM(01737, 02656);
	SaveEMEM(01740, 014434);
	SaveEMEM(01741, 067336);
	SaveEMEM(01742, 043072);
	SaveEMEM(01743, 015154);
	SaveEMEM(01744, 015512);
	SaveEMEM(01745, 063006);
	SaveEMEM(01746, 064677);
	SaveEMEM(01747, 06242);
	SaveEMEM(01750, 0131);
	SaveEMEM(01751, 06706);
	SaveEMEM(01752, 020152);

	//RANGEVAR
	DoubleToBuffer(1.111111111e-5, -12, iTemp, iTemp2);
	SaveEMEM(01766, iTemp);
	SaveEMEM(01767, iTemp2);
	//RATEVAR
	DoubleToBuffer(1.877777778e-5, -12, iTemp, iTemp2);
	SaveEMEM(01770, iTemp);
	SaveEMEM(01771, iTemp2);
	//RVARMIN
	iTemp = SingleToBuffer(66.0, 12);
	SaveEMEM(01772, iTemp);
	//VVARMIN
	iTemp = SingleToBuffer(1.7445e-6, -12);
	SaveEMEM(01773, iTemp);

	//WSURFPOS
	SaveEMEM(02006, 0);
	//WSURFVEL
	SaveEMEM(02007, 0);

	//SHAFTVAR
	iTemp = SingleToBuffer(BLOCKII.SHAFTVAR*1e-6, -12);
	SaveEMEM(02010, iTemp);
	//TRUNVAR
	iTemp = SingleToBuffer(BLOCKII.TRUNVAR*1e-6, -12);
	SaveEMEM(02011, iTemp);

	//RLS
	RLS = r_from_latlong(LSLat*RAD, LSLng*RAD, LSAlt, BODY_MOON, 1);
	DoubleToBuffer(RLS.x, 27, iTemp, iTemp2);
	SaveEMEM(02020, iTemp);
	SaveEMEM(02021, iTemp2);
	DoubleToBuffer(RLS.y, 27, iTemp, iTemp2);
	SaveEMEM(02022, iTemp);
	SaveEMEM(02023, iTemp2);
	DoubleToBuffer(RLS.z, 27, iTemp, iTemp2);
	SaveEMEM(02024, iTemp);
	SaveEMEM(02025, iTemp2);

	//TLAND
	DoubleToBuffer(BLOCKII.TLAND*3600.0*100.0, 28, iTemp, iTemp2);
	SaveEMEM(02026, iTemp);
	SaveEMEM(02027, iTemp2);

	//VELBIAS
	DoubleToBuffer(2.5*0.3048 / 100.0, 6, iTemp, iTemp2);
	SaveEMEM(02400, iTemp);
	SaveEMEM(02401, iTemp2);

	BLOCKII.RBRFG.x = -1773.725;
	BLOCKII.RAPFG.x = 94.91910105;
	BLOCKII.RBRFG.z = -14488.0269;
	BLOCKII.RAPFG.z = -15.72079987;
	BLOCKII.VBRFG.x = -168.1064567;
	BLOCKII.VAPFG.x = 2.083579987;
	BLOCKII.VBRFG.z = -77.6143668;
	BLOCKII.VAPFG.z = 0.8303187992;
	BLOCKII.ABRFG.x = -0.6472360236;
	BLOCKII.AAPFG.x = 0.5402850066;
	BLOCKII.ABRFG.z = -8.41438189;
	BLOCKII.AAPFG.z = -0.2354229987;
	BLOCKII.VBRFG_star = -1397.058596;
	BLOCKII.VAPFG_star = 14.94573786;
	BLOCKII.ABRFG_star = -50.48629265;
	BLOCKII.AAPFG_star = -1.412537992;
	BLOCKII.JBRFG_star = 0.008257294947;
	BLOCKII.JAPFG_star = 0.04509242126;

	BLOCKII.TCGFBRAK = 30.0;
	BLOCKII.TCGIBRAK = 900.0;
	BLOCKII.TCGFAPPR = 6.0;
	BLOCKII.TCGIAPPR = 200.0;

	BLOCKII.DELQFIX = 100.0;

	//RBRFGX
	DoubleToBuffer(BLOCKII.RBRFG.x*FT2M, 24, iTemp, iTemp2);
	SaveEMEM(02402, iTemp);
	SaveEMEM(02403, iTemp2);
	//RAPFGX
	DoubleToBuffer(BLOCKII.RAPFG.x*FT2M, 24, iTemp, iTemp2);
	SaveEMEM(02404, iTemp);
	SaveEMEM(02405, iTemp2);
	//RBRFGZ
	DoubleToBuffer(BLOCKII.RBRFG.z*FT2M, 24, iTemp, iTemp2);
	SaveEMEM(02406, iTemp);
	SaveEMEM(02407, iTemp2);
	//RAPFGZ
	DoubleToBuffer(BLOCKII.RAPFG.z*FT2M, 24, iTemp, iTemp2);
	SaveEMEM(02410, iTemp);
	SaveEMEM(02411, iTemp2);
	//VBRFGX
	DoubleToBuffer(BLOCKII.VBRFG.x*FT2M / 100.0, 10, iTemp, iTemp2);
	SaveEMEM(02412, iTemp);
	SaveEMEM(02413, iTemp2);
	//VAPFGX
	DoubleToBuffer(BLOCKII.VAPFG.x*FT2M / 100.0, 10, iTemp, iTemp2);
	SaveEMEM(02414, iTemp);
	SaveEMEM(02415, iTemp2);
	//VBRFGZ
	DoubleToBuffer(BLOCKII.VBRFG.z*FT2M / 100.0, 10, iTemp, iTemp2);
	SaveEMEM(02416, iTemp);
	SaveEMEM(02417, iTemp2);
	//VAPFGZ
	DoubleToBuffer(BLOCKII.VAPFG.z*FT2M / 100.0, 10, iTemp, iTemp2);
	SaveEMEM(02420, iTemp);
	SaveEMEM(02421, iTemp2);
	//ABRFGX
	DoubleToBuffer(BLOCKII.ABRFG.x*FT2M / 10000.0, -4, iTemp, iTemp2);
	SaveEMEM(02422, iTemp);
	SaveEMEM(02423, iTemp2);
	//AAPFGX
	DoubleToBuffer(BLOCKII.AAPFG.x*FT2M / 10000.0, -4, iTemp, iTemp2);
	SaveEMEM(02424, iTemp);
	SaveEMEM(02425, iTemp2);
	//ABRFGZ
	DoubleToBuffer(BLOCKII.ABRFG.z*FT2M / 10000.0, -4, iTemp, iTemp2);
	SaveEMEM(02426, iTemp);
	SaveEMEM(02427, iTemp2);
	//AAPFGZ
	DoubleToBuffer(BLOCKII.AAPFG.z*FT2M / 10000.0, -4, iTemp, iTemp2);
	SaveEMEM(02430, iTemp);
	SaveEMEM(02431, iTemp2);
	//VBRFG*
	DoubleToBuffer(BLOCKII.VBRFG_star*FT2M / 100.0, 13, iTemp, iTemp2);
	SaveEMEM(02432, iTemp);
	SaveEMEM(02433, iTemp2);
	//VAPFG*
	DoubleToBuffer(BLOCKII.VAPFG_star*FT2M / 100.0, 13, iTemp, iTemp2);
	SaveEMEM(02434, iTemp);
	SaveEMEM(02435, iTemp2);
	//ABRFG*
	DoubleToBuffer(BLOCKII.ABRFG_star*FT2M / 10000.0, -4, iTemp, iTemp2);
	SaveEMEM(02436, iTemp);
	SaveEMEM(02437, iTemp2);
	//AAPFG*
	DoubleToBuffer(BLOCKII.AAPFG_star*FT2M / 10000.0, -4, iTemp, iTemp2);
	SaveEMEM(02440, iTemp);
	SaveEMEM(02441, iTemp2);
	//JBRFG*
	DoubleToBuffer(BLOCKII.JBRFG_star*FT2M / 1000000.0, -21, iTemp, iTemp2);
	SaveEMEM(02442, iTemp);
	SaveEMEM(02443, iTemp2);
	//JAPFG*
	DoubleToBuffer(BLOCKII.JAPFG_star*FT2M / 1000000.0, -21, iTemp, iTemp2);
	SaveEMEM(02444, iTemp);
	SaveEMEM(02445, iTemp2);
	//GAINBRAK
	SaveEMEM(02446, 037777);
	SaveEMEM(02447, 037777);
	//GAINAPPR
	SaveEMEM(02450, 0);
	SaveEMEM(02451, 0);
	//TCGFBRAK
	iTemp = SingleToBuffer(BLOCKII.TCGFBRAK*100.0, 17);
	SaveEMEM(02452, iTemp);
	//TCGIBRAK
	iTemp = SingleToBuffer(BLOCKII.TCGIBRAK*100.0, 17);
	SaveEMEM(02453, iTemp);
	//TCGFAPPR
	iTemp = SingleToBuffer(BLOCKII.TCGFAPPR*100.0, 17);
	SaveEMEM(02454, iTemp);
	//TCGIAPPR
	iTemp = SingleToBuffer(BLOCKII.TCGIAPPR*100.0, 17);
	SaveEMEM(02455, iTemp);

	//VIGN
	DoubleToBuffer(BLOCKII.VIGN*FT2M / 100.0, 10, iTemp, iTemp2);
	SaveEMEM(02456, iTemp);
	SaveEMEM(02457, iTemp2);
	//RIGNX
	DoubleToBuffer(BLOCKII.RIGNX*0.3048, 24, iTemp, iTemp2);
	SaveEMEM(02460, iTemp);
	SaveEMEM(02461, iTemp2);
	//RIGNZ
	DoubleToBuffer(BLOCKII.RIGNZ*0.3048, 24, iTemp, iTemp2);
	SaveEMEM(02462, iTemp);
	SaveEMEM(02463, iTemp2);
	//KIGNX/B4
	DoubleToBuffer(BLOCKII.KIGNXB4, 4, iTemp, iTemp2);
	SaveEMEM(02464, iTemp);
	SaveEMEM(02465, iTemp2);
	//KIGNY/B8
	DoubleToBuffer(BLOCKII.KIGNYB8 / 0.3048, -16, iTemp, iTemp2);
	SaveEMEM(02466, iTemp);
	SaveEMEM(02467, iTemp2);
	//KIGNV/B4
	DoubleToBuffer(BLOCKII.KIGNVB4*100.0, 18, iTemp, iTemp2);
	SaveEMEM(02470, iTemp);
	SaveEMEM(02471, iTemp2);
	//LOWCRIT
	SaveEMEM(02472, 04114);
	//HIGHCRIT
	SaveEMEM(02473, 04454);

	//TAUHZ
	iTemp = SingleToBuffer(5.0*100.0, 10);
	SaveEMEM(02474, iTemp);
	//QHZ
	SaveEMEM(02475, 031463);
	//AHZLIM
	SaveEMEM(02476, 037);
	//TEXTRA (Time to achieve P66HZ command)
	iTemp = SingleToBuffer(0.75*100.0, 9);
	SaveEMEM(02477, iTemp);

	//DELQFIX
	DoubleToBuffer(BLOCKII.DELQFIX*FT2M, 24, iTemp, iTemp2);
	SaveEMEM(02500, iTemp);
	SaveEMEM(02501, iTemp2);
	//LRVMAX
	iTemp = SingleToBuffer(2500.0*0.3048 / 100.0, 7);
	SaveEMEM(02502, iTemp);
	//LRVF
	iTemp = SingleToBuffer(200.0*0.3048 / 100.0, 7);
	SaveEMEM(02503, iTemp);
	//LRWVZ
	iTemp = SingleToBuffer(0.3, 0);
	SaveEMEM(02504, iTemp);
	//LRWVY
	iTemp = SingleToBuffer(0.3, 0);
	SaveEMEM(02505, iTemp);
	//LRWVX
	iTemp = SingleToBuffer(0.3, 0);
	SaveEMEM(02506, iTemp);
	//LRWVFZ
	iTemp = SingleToBuffer(0.2, 0);
	SaveEMEM(02507, iTemp);
	//LRWVFY
	iTemp = SingleToBuffer(0.2, 0);
	SaveEMEM(02510, iTemp);
	//LRWVFX
	iTemp = SingleToBuffer(0.2, 0);
	SaveEMEM(02511, iTemp);
	//LRWVFF
	iTemp = SingleToBuffer(0.1, 0);
	SaveEMEM(02512, iTemp);

	//ABSC0
	iTemp = SingleToBuffer(BLOCKII.ABSC[0] * FT2M, 18);
	SaveEMEM(02513, iTemp);
	//ABSC1
	iTemp = SingleToBuffer(BLOCKII.ABSC[1] * FT2M, 18);
	SaveEMEM(02514, iTemp);
	//ABSC2
	iTemp = SingleToBuffer(BLOCKII.ABSC[2] * FT2M, 18);
	SaveEMEM(02515, iTemp);
	//ABSC3
	iTemp = SingleToBuffer(BLOCKII.ABSC[3] * FT2M, 18);
	SaveEMEM(02516, iTemp);
	//ABSC4
	iTemp = SingleToBuffer(BLOCKII.ABSC[4] * FT2M, 18);
	SaveEMEM(02517, iTemp);
	//SLOPE0
	iTemp = SingleToBuffer(BLOCKII.SLOPE[0], 6);
	SaveEMEM(02520, iTemp);
	//SLOPE1
	iTemp = SingleToBuffer(BLOCKII.SLOPE[1], 6);
	SaveEMEM(02521, iTemp);
	//SLOPE2
	iTemp = SingleToBuffer(BLOCKII.SLOPE[2], 6);
	SaveEMEM(02522, iTemp);
	//SLOPE3
	iTemp = SingleToBuffer(BLOCKII.SLOPE[3], 6);
	SaveEMEM(02523, iTemp);
	//SLOPE4
	iTemp = SingleToBuffer(BLOCKII.SLOPE[4], 6);
	SaveEMEM(02524, iTemp);

	//RODSCALE
	iTemp = SingleToBuffer(1.0*0.3048 / 100.0, -3);
	SaveEMEM(02525, iTemp);
	//ROHZSCAL
	iTemp = SingleToBuffer(1.0*0.3048 / 100.0, -3);
	SaveEMEM(02526, iTemp);

	//TAURODL (Little Tau for P66ROD outside DB)
	iTemp = SingleToBuffer(1.5*100.0, 14);
	SaveEMEM(02527, iTemp);
	//TAURODB (Bigger Tau for P66ROD inside DB)
	iTemp = SingleToBuffer(3.0*100.0, 14);
	SaveEMEM(02530, iTemp);
	//VERCRIT (Velocity error criterion for P66ROD)
	SaveEMEM(02531, 037777);

	//MINFORCE
	SaveEMEM(02532, 0534);
	//MAXFORCE
	SaveEMEM(02533, 04274);

	//J1PARM
	DoubleToBuffer(BLOCKII.J1PARM*FT2M, 23, iTemp, iTemp2);
	SaveEMEM(02534, iTemp);
	SaveEMEM(02535, iTemp2);
	//K1PARM
	DoubleToBuffer(BLOCKII.K1PARM*FT2M*PI2, 23, iTemp, iTemp2);
	SaveEMEM(02536, iTemp);
	SaveEMEM(02537, iTemp2);
	//J2PARM
	DoubleToBuffer(BLOCKII.J2PARM*FT2M, 23, iTemp, iTemp2);
	SaveEMEM(02540, iTemp);
	SaveEMEM(02541, iTemp2);
	//K2PARM
	DoubleToBuffer(BLOCKII.K2PARM*FT2M*PI2, 23, iTemp, iTemp2);
	SaveEMEM(02542, iTemp);
	SaveEMEM(02543, iTemp2);
	//THETCRIT
	DoubleToBuffer(BLOCKII.THETCRIT / 360.0, 0, iTemp, iTemp2);
	SaveEMEM(02544, iTemp);
	SaveEMEM(02545, iTemp2);
	//RAMIN
	DoubleToBuffer(BLOCKII.RAMIN*FT2M, 24, iTemp, iTemp2);
	SaveEMEM(02546, iTemp);
	SaveEMEM(02547, iTemp2);
	//YLIM
	DoubleToBuffer(8.2*1852.0, 24, iTemp, iTemp2);
	SaveEMEM(02550, iTemp);
	SaveEMEM(02551, iTemp2);
	//ABTRDOT
	DoubleToBuffer(19.5*0.3048 / 100.0, 7, iTemp, iTemp2);
	SaveEMEM(02552, iTemp);
	SaveEMEM(02553, iTemp2);
	//COSTHET1
	SaveEMEM(02554, 0);
	SaveEMEM(02555, 0);
	//COSTHET2
	SaveEMEM(02556, 06733);
	SaveEMEM(02557, 07535);

	//DLAND
	SaveEMEM(02612, 0);
	SaveEMEM(02613, 0);
	SaveEMEM(02614, 0);
	SaveEMEM(02615, 0);
	SaveEMEM(02616, 0);
	SaveEMEM(02617, 0);

	//IGNAOSQ
	iTemp = SingleToBuffer(BLOCKII.IGNAOSQ / 360.0, -2);
	SaveEMEM(03012, iTemp);
	//IGNAOSR
	iTemp = SingleToBuffer(BLOCKII.IGNAOSR / 360.0, -2);
	SaveEMEM(03013, iTemp);

	//DOWNTORK
	SaveEMEM(03113, 0);
	SaveEMEM(03114, 0);
	SaveEMEM(03115, 0);
	SaveEMEM(03116, 0);
	SaveEMEM(03117, 0);
	SaveEMEM(03120, 0);

	//AGSK
	DoubleToBuffer(BLOCKII.AGSK*3600.0*100.0, 28, iTemp, iTemp2);
	SaveEMEM(03371, iTemp);
	SaveEMEM(03372, iTemp2);

	//AZBIAS
	iTemp = SingleToBuffer(0.0 / 360.0, -1);
	SaveEMEM(03373, iTemp);

	//AOTAZ
	iTemp = SingleToBuffer(-60.0 / 360.0, -1, true);
	SaveEMEM(03404, iTemp);
	iTemp = SingleToBuffer(0.0 / 360.0, -1, true);
	SaveEMEM(03405, iTemp);
	iTemp = SingleToBuffer(60.0 / 360.0, -1, true);
	SaveEMEM(03406, iTemp);
	iTemp = SingleToBuffer(120.0 / 360.0, -1, true);
	SaveEMEM(03407, iTemp);
	iTemp = SingleToBuffer(-180.0 / 360.0, -1, true);
	SaveEMEM(03410, iTemp);
	iTemp = SingleToBuffer(-120.0 / 360.0, -1, true);
	SaveEMEM(03411, iTemp);
	//AOTEL
	iTemp = SingleToBuffer(45.0 / 360.0, -1);
	SaveEMEM(03412, iTemp);
	iTemp = SingleToBuffer(45.0 / 360.0, -1);
	SaveEMEM(03413, iTemp);
	iTemp = SingleToBuffer(45.0 / 360.0, -1);
	SaveEMEM(03414, iTemp);
	iTemp = SingleToBuffer(45.0 / 360.0, -1);
	SaveEMEM(03415, iTemp);
	iTemp = SingleToBuffer(45.0 / 360.0, -1);
	SaveEMEM(03416, iTemp);
	iTemp = SingleToBuffer(45.0 / 360.0, -1);
	SaveEMEM(03417, iTemp);

	//LRHMAX
	iTemp = SingleToBuffer(50000.0*0.3048, 14);
	SaveEMEM(03420, iTemp);
	//LRWH
	iTemp = SingleToBuffer(0.35, 0);
	SaveEMEM(03421, iTemp);
	//ZOOMTIME
	iTemp = SingleToBuffer(26.0*100.0, 14);
	SaveEMEM(03422, iTemp);
	//TENDBRAK
	iTemp = SingleToBuffer(62.0*100.0, 17);
	SaveEMEM(03423, iTemp);
	//TENDAPPR
	iTemp = SingleToBuffer(12.0*100.0, 17);
	SaveEMEM(03424, iTemp);
	//DELTTFAP
	iTemp = SingleToBuffer(BLOCKII.DELTTFAP*100.0, 17);
	SaveEMEM(03425, iTemp);
	//LEADTIME
	iTemp = SingleToBuffer(-2.2*100.0, 17);
	SaveEMEM(03426, iTemp);
	//RPCRTIME
	iTemp = SingleToBuffer(62.0*100.0, 17);
	SaveEMEM(03427, iTemp);
	//RPCRTQSW
	iTemp = SingleToBuffer(-1.0, 1);
	SaveEMEM(03430, iTemp);
	//TNEWA
	SaveEMEM(03431, 020000);
	SaveEMEM(03432, 0);
}

void AGCPadloadGenerator::Luminary210Padload()
{
	LGCDefaults(false);

	//CHANBKUP
	SaveEMEM(0374, 00011);

	//LRWH1
	iTemp = SingleToBuffer(0.35, 0);
	SaveEMEM(01315, iTemp);

	//REFSMMAT
	SaveEMEM(01731, 011244);
	SaveEMEM(01732, 012736);
	SaveEMEM(01733, 07474);
	SaveEMEM(01734, 035041);
	SaveEMEM(01735, 012423);
	SaveEMEM(01736, 037124);
	SaveEMEM(01737, 073772);
	SaveEMEM(01740, 064052);
	SaveEMEM(01741, 065262);
	SaveEMEM(01742, 077044);
	SaveEMEM(01743, 013176);
	SaveEMEM(01744, 06616);
	SaveEMEM(01745, 014275);
	SaveEMEM(01746, 032203);
	SaveEMEM(01747, 066634);
	SaveEMEM(01750, 076474);
	SaveEMEM(01751, 073551);
	SaveEMEM(01752, 041767);

	BLOCKII.RBRFG.x = -3118.358793;
	BLOCKII.RAPFG.x = 158.5;
	BLOCKII.RBRFG.z = -11741.44094;
	BLOCKII.RAPFG.z = -27.3554;
	BLOCKII.VBRFG.x = -196.46916;
	BLOCKII.VAPFG.x = -3.534759843;
	BLOCKII.VBRFG.z = -166.7599672;
	BLOCKII.VAPFG.z = 0.0249505;
	BLOCKII.ABRFG.x = -0.7182481299;
	BLOCKII.AAPFG.x = 0.07717830052;
	BLOCKII.ABRFG.z = -8.302450459;
	BLOCKII.AAPFG.z = -0.5896270013;
	BLOCKII.VBRFG_star = -3001.679987;
	BLOCKII.VAPFG_star = 0.4491089895;
	BLOCKII.ABRFG_star = -49.81470144;
	BLOCKII.AAPFG_star = -3.537762139;
	BLOCKII.JBRFG_star = -0.01512365879;
	BLOCKII.JAPFG_star = 0.04317359908;

	BLOCKII.TCGFBRAK = 30.0;
	BLOCKII.TCGIBRAK = 900.0;
	BLOCKII.TCGFAPPR = 6.0;
	BLOCKII.TCGIAPPR = 200.0;

	BLOCKII.DELQFIX = 100.0;

	DescentConstants14_17();

	//N26/PRI
	SaveEMEM(02371, 0);
	//N26/2CAD
	SaveEMEM(02372, 0);
	SaveEMEM(02373, 0);
}

void AGCPadloadGenerator::DescentConstants14_17()
{
	//Same addresses from Apollo 14 to 17

	//FLAGWRD3
	SaveEMEM(077, 012000);
	//FLAGWRD8
	SaveEMEM(0104, 06000);
	//FLAGWRD10
	SaveEMEM(0106, 0);

	//MASS
	DoubleToBuffer(BLOCKII.TotalMass*LBS2KG, 16, iTemp, iTemp2);
	SaveEMEM(01243, iTemp);
	SaveEMEM(01244, 0);

	//LEMMASS
	iTemp = SingleToBuffer(BLOCKII.LMMass*LBS2KG, 16);
	SaveEMEM(01326, iTemp);
	//CSMMASS
	iTemp = SingleToBuffer(BLOCKII.CSMMass*LBS2KG, 16);
	SaveEMEM(01327, iTemp);

	R2Model(01347);

	//ELBIAS
	SaveEMEM(01353, 0);

	//TOOFEW
	SaveEMEM(01354, 03);

	//TETCSM
	SaveEMEM(01570, 037777);
	SaveEMEM(01571, 037777);
	//TETLEM
	SaveEMEM(01642, 037777);
	SaveEMEM(01643, 037777);

	//RANGEVAR
	DoubleToBuffer(1.111111111e-5, -12, iTemp, iTemp2);
	SaveEMEM(01766, iTemp);
	SaveEMEM(01767, iTemp2);
	//RATEVAR
	DoubleToBuffer(1.877777778e-5, -12, iTemp, iTemp2);
	SaveEMEM(01770, iTemp);
	SaveEMEM(01771, iTemp2);
	//RVARMIN
	iTemp = SingleToBuffer(66.0, 12);
	SaveEMEM(01772, iTemp);
	//VVARMIN
	iTemp = SingleToBuffer(1.7445e-6, -12);
	SaveEMEM(01773, iTemp);

	//WSURFPOS
	SaveEMEM(02006, 0);
	//WSURFVEL
	SaveEMEM(02007, 0);

	//SHAFTVAR
	iTemp = SingleToBuffer(BLOCKII.SHAFTVAR*1e-6, -12);
	SaveEMEM(02010, iTemp);
	//TRUNVAR
	iTemp = SingleToBuffer(BLOCKII.TRUNVAR*1e-6, -12);
	SaveEMEM(02011, iTemp);

	//RLS
	RLS = r_from_latlong(LSLat*RAD, LSLng*RAD, LSAlt, BODY_MOON, 1);
	DoubleToBuffer(RLS.x, 27, iTemp, iTemp2);
	SaveEMEM(02020, iTemp);
	SaveEMEM(02021, iTemp2);
	DoubleToBuffer(RLS.y, 27, iTemp, iTemp2);
	SaveEMEM(02022, iTemp);
	SaveEMEM(02023, iTemp2);
	DoubleToBuffer(RLS.z, 27, iTemp, iTemp2);
	SaveEMEM(02024, iTemp);
	SaveEMEM(02025, iTemp2);

	//TLAND
	DoubleToBuffer(BLOCKII.TLAND*3600.0*100.0, 28, iTemp, iTemp2);
	SaveEMEM(02026, iTemp);
	SaveEMEM(02027, iTemp2);

	//VELBIAS
	DoubleToBuffer(2.5*0.3048 / 100.0, 6, iTemp, iTemp2);
	SaveEMEM(02400, iTemp);
	SaveEMEM(02401, iTemp2);

	//RBRFGX
	DoubleToBuffer(BLOCKII.RBRFG.x*FT2M, 24, iTemp, iTemp2);
	SaveEMEM(02402, iTemp);
	SaveEMEM(02403, iTemp2);
	//RAPFGX
	DoubleToBuffer(BLOCKII.RAPFG.x*FT2M, 24, iTemp, iTemp2);
	SaveEMEM(02404, iTemp);
	SaveEMEM(02405, iTemp2);
	//RBRFGZ
	DoubleToBuffer(BLOCKII.RBRFG.z*FT2M, 24, iTemp, iTemp2);
	SaveEMEM(02406, iTemp);
	SaveEMEM(02407, iTemp2);
	//RAPFGZ
	DoubleToBuffer(BLOCKII.RAPFG.z*FT2M, 24, iTemp, iTemp2);
	SaveEMEM(02410, iTemp);
	SaveEMEM(02411, iTemp2);
	//VBRFGX
	DoubleToBuffer(BLOCKII.VBRFG.x*FT2M / 100.0, 10, iTemp, iTemp2);
	SaveEMEM(02412, iTemp);
	SaveEMEM(02413, iTemp2);
	//VAPFGX
	DoubleToBuffer(BLOCKII.VAPFG.x*FT2M / 100.0, 10, iTemp, iTemp2);
	SaveEMEM(02414, iTemp);
	SaveEMEM(02415, iTemp2);
	//VBRFGZ
	DoubleToBuffer(BLOCKII.VBRFG.z*FT2M / 100.0, 10, iTemp, iTemp2);
	SaveEMEM(02416, iTemp);
	SaveEMEM(02417, iTemp2);
	//VAPFGZ
	DoubleToBuffer(BLOCKII.VAPFG.z*FT2M / 100.0, 10, iTemp, iTemp2);
	SaveEMEM(02420, iTemp);
	SaveEMEM(02421, iTemp2);
	//ABRFGX
	DoubleToBuffer(BLOCKII.ABRFG.x*FT2M / 10000.0, -4, iTemp, iTemp2);
	SaveEMEM(02422, iTemp);
	SaveEMEM(02423, iTemp2);
	//AAPFGX
	DoubleToBuffer(BLOCKII.AAPFG.x*FT2M / 10000.0, -4, iTemp, iTemp2);
	SaveEMEM(02424, iTemp);
	SaveEMEM(02425, iTemp2);
	//ABRFGZ
	DoubleToBuffer(BLOCKII.ABRFG.z*FT2M / 10000.0, -4, iTemp, iTemp2);
	SaveEMEM(02426, iTemp);
	SaveEMEM(02427, iTemp2);
	//AAPFGZ
	DoubleToBuffer(BLOCKII.AAPFG.z*FT2M / 10000.0, -4, iTemp, iTemp2);
	SaveEMEM(02430, iTemp);
	SaveEMEM(02431, iTemp2);
	//VBRFG*
	DoubleToBuffer(BLOCKII.VBRFG_star*FT2M / 100.0, 13, iTemp, iTemp2);
	SaveEMEM(02432, iTemp);
	SaveEMEM(02433, iTemp2);
	//VAPFG*
	DoubleToBuffer(BLOCKII.VAPFG_star*FT2M / 100.0, 13, iTemp, iTemp2);
	SaveEMEM(02434, iTemp);
	SaveEMEM(02435, iTemp2);
	//ABRFG*
	DoubleToBuffer(BLOCKII.ABRFG_star*FT2M / 10000.0, -4, iTemp, iTemp2);
	SaveEMEM(02436, iTemp);
	SaveEMEM(02437, iTemp2);
	//AAPFG*
	DoubleToBuffer(BLOCKII.AAPFG_star*FT2M / 10000.0, -4, iTemp, iTemp2);
	SaveEMEM(02440, iTemp);
	SaveEMEM(02441, iTemp2);
	//JBRFG*
	DoubleToBuffer(BLOCKII.JBRFG_star*FT2M / 1000000.0, -21, iTemp, iTemp2);
	SaveEMEM(02442, iTemp);
	SaveEMEM(02443, iTemp2);
	//JAPFG*
	DoubleToBuffer(BLOCKII.JAPFG_star*FT2M / 1000000.0, -21, iTemp, iTemp2);
	SaveEMEM(02444, iTemp);
	SaveEMEM(02445, iTemp2);
	//GAINBRAK
	SaveEMEM(02446, 037777);
	SaveEMEM(02447, 037777);
	//GAINAPPR
	SaveEMEM(02450, 0);
	SaveEMEM(02451, 0);
	//TCGFBRAK
	iTemp = SingleToBuffer(BLOCKII.TCGFBRAK*100.0, 17);
	SaveEMEM(02452, iTemp);
	//TCGIBRAK
	iTemp = SingleToBuffer(BLOCKII.TCGIBRAK*100.0, 17);
	SaveEMEM(02453, iTemp);
	//TCGFAPPR
	iTemp = SingleToBuffer(BLOCKII.TCGFAPPR*100.0, 17);
	SaveEMEM(02454, iTemp);
	//TCGIAPPR
	iTemp = SingleToBuffer(BLOCKII.TCGIAPPR*100.0, 17);
	SaveEMEM(02455, iTemp);

	//VIGN
	DoubleToBuffer(BLOCKII.VIGN*FT2M / 100.0, 10, iTemp, iTemp2);
	SaveEMEM(02456, iTemp);
	SaveEMEM(02457, iTemp2);
	//RIGNX
	DoubleToBuffer(BLOCKII.RIGNX*0.3048, 24, iTemp, iTemp2);
	SaveEMEM(02460, iTemp);
	SaveEMEM(02461, iTemp2);
	//RIGNZ
	DoubleToBuffer(BLOCKII.RIGNZ*0.3048, 24, iTemp, iTemp2);
	SaveEMEM(02462, iTemp);
	SaveEMEM(02463, iTemp2);
	//KIGNX/B4
	DoubleToBuffer(BLOCKII.KIGNXB4, 4, iTemp, iTemp2);
	SaveEMEM(02464, iTemp);
	SaveEMEM(02465, iTemp2);
	//KIGNY/B8
	DoubleToBuffer(BLOCKII.KIGNYB8 / 0.3048, -16, iTemp, iTemp2);
	SaveEMEM(02466, iTemp);
	SaveEMEM(02467, iTemp2);
	//KIGNV/B4
	DoubleToBuffer(BLOCKII.KIGNVB4*100.0, 18, iTemp, iTemp2);
	SaveEMEM(02470, iTemp);
	SaveEMEM(02471, iTemp2);
	//LOWCRIT
	SaveEMEM(02472, 04114);
	//HIGHCRIT
	SaveEMEM(02473, 04454);

	//TAUHZ
	iTemp = SingleToBuffer(5.0*100.0, 11);
	SaveEMEM(02474, iTemp);
	//QHZ
	SaveEMEM(02475, 014632);
	//AHZLIM
	SaveEMEM(02476, 017);
	//2LATE466
	DoubleToBuffer(1.50*100.0, 28, iTemp, iTemp2);
	SaveEMEM(02477, iTemp);
	SaveEMEM(02500, iTemp2);
	//DELQFIX
	DoubleToBuffer(BLOCKII.DELQFIX*FT2M, 24, iTemp, iTemp2);
	SaveEMEM(02503, iTemp);
	SaveEMEM(02504, iTemp2);
	//LRVMAX
	iTemp = SingleToBuffer(2500.0*0.3048 / 100.0, 7);
	SaveEMEM(02511, iTemp);
	//LRVF
	iTemp = SingleToBuffer(200.0*0.3048 / 100.0, 7);
	SaveEMEM(02512, iTemp);
	//LRWVZ
	iTemp = SingleToBuffer(0.3, 0);
	SaveEMEM(02513, iTemp);
	//LRWVY
	iTemp = SingleToBuffer(0.3, 0);
	SaveEMEM(02514, iTemp);
	//LRWVX
	iTemp = SingleToBuffer(0.3, 0);
	SaveEMEM(02515, iTemp);
	//LRWVFZ
	iTemp = SingleToBuffer(0.2, 0);
	SaveEMEM(02516, iTemp);
	//LRWVFY
	iTemp = SingleToBuffer(0.2, 0);
	SaveEMEM(02517, iTemp);
	//LRWVFX
	iTemp = SingleToBuffer(0.2, 0);
	SaveEMEM(02520, iTemp);
	//LRWVFF
	iTemp = SingleToBuffer(0.1, 0);
	SaveEMEM(02521, iTemp);

	//ABSC0
	iTemp = SingleToBuffer(BLOCKII.ABSC[0] * FT2M, 18);
	SaveEMEM(02522, iTemp);
	//ABSC1
	iTemp = SingleToBuffer(BLOCKII.ABSC[1] * FT2M, 18);
	SaveEMEM(02523, iTemp);
	//ABSC2
	iTemp = SingleToBuffer(BLOCKII.ABSC[2] * FT2M, 18);
	SaveEMEM(02524, iTemp);
	//ABSC3
	iTemp = SingleToBuffer(BLOCKII.ABSC[3] * FT2M, 18);
	SaveEMEM(02525, iTemp);
	//ABSC4
	iTemp = SingleToBuffer(BLOCKII.ABSC[4] * FT2M, 18);
	SaveEMEM(02526, iTemp);
	//SLOPE0
	iTemp = SingleToBuffer(BLOCKII.SLOPE[0], 6);
	SaveEMEM(02527, iTemp);
	//SLOPE1
	iTemp = SingleToBuffer(BLOCKII.SLOPE[1], 6);
	SaveEMEM(02530, iTemp);
	//SLOPE2
	iTemp = SingleToBuffer(BLOCKII.SLOPE[2], 6);
	SaveEMEM(02531, iTemp);
	//SLOPE3
	iTemp = SingleToBuffer(BLOCKII.SLOPE[3], 6);
	SaveEMEM(02532, iTemp);
	//SLOPE4
	iTemp = SingleToBuffer(BLOCKII.SLOPE[4], 6);
	SaveEMEM(02533, iTemp);

	//RODSCALE
	iTemp = SingleToBuffer(1.0*0.3048 / 100.0, -7);
	SaveEMEM(02534, iTemp);
	//TAUROD
	DoubleToBuffer(1.5*100.0, 9, iTemp, iTemp2);
	SaveEMEM(02535, iTemp);
	SaveEMEM(02536, iTemp2);
	//LAG/TAU
	DoubleToBuffer(0.23333, 0, iTemp, iTemp2);
	SaveEMEM(02537, iTemp);
	SaveEMEM(02540, iTemp2);
	//MINFORCE
	SaveEMEM(02541, 01);
	SaveEMEM(02542, 027631);
	//MAXFORCE
	SaveEMEM(02543, 013);
	SaveEMEM(02544, 06551);

	//J1PARM
	DoubleToBuffer(BLOCKII.J1PARM*FT2M, 23, iTemp, iTemp2);
	SaveEMEM(02545, iTemp);
	SaveEMEM(02546, iTemp2);
	//K1PARM
	DoubleToBuffer(BLOCKII.K1PARM*FT2M*PI2, 23, iTemp, iTemp2);
	SaveEMEM(02547, iTemp);
	SaveEMEM(02550, iTemp2);
	//J2PARM
	DoubleToBuffer(BLOCKII.J2PARM*FT2M, 23, iTemp, iTemp2);
	SaveEMEM(02551, iTemp);
	SaveEMEM(02552, iTemp2);
	//K2PARM
	DoubleToBuffer(BLOCKII.K2PARM*FT2M*PI2, 23, iTemp, iTemp2);
	SaveEMEM(02553, iTemp);
	SaveEMEM(02554, iTemp2);
	//THETCRIT
	DoubleToBuffer(BLOCKII.THETCRIT / 360.0, 0, iTemp, iTemp2);
	SaveEMEM(02555, iTemp);
	SaveEMEM(02556, iTemp2);
	//RAMIN
	DoubleToBuffer(BLOCKII.RAMIN*FT2M, 24, iTemp, iTemp2);
	SaveEMEM(02557, iTemp);
	SaveEMEM(02560, iTemp2);
	//YLIM
	DoubleToBuffer(8.2*1852.0, 24, iTemp, iTemp2);
	SaveEMEM(02561, iTemp);
	SaveEMEM(02562, iTemp2);
	//ABTRDOT
	DoubleToBuffer(19.5*0.3048 / 100.0, 7, iTemp, iTemp2);
	SaveEMEM(02563, iTemp);
	SaveEMEM(02564, iTemp2);
	//COSTHET1
	SaveEMEM(02565, 0);
	SaveEMEM(02566, 0);
	//COSTHET2
	SaveEMEM(02567, 06733);
	SaveEMEM(02570, 07535);

	//DLAND
	SaveEMEM(02631, 0);
	SaveEMEM(02632, 0);
	SaveEMEM(02633, 0);
	SaveEMEM(02634, 0);
	SaveEMEM(02635, 0);
	SaveEMEM(02636, 0);

	//IGNAOSQ
	iTemp = SingleToBuffer(BLOCKII.IGNAOSQ / 360.0, -2);
	SaveEMEM(03012, iTemp);
	//IGNAOSR
	iTemp = SingleToBuffer(BLOCKII.IGNAOSR / 360.0, -2);
	SaveEMEM(03013, iTemp);

	//DOWNTORK
	SaveEMEM(03113, 0);
	SaveEMEM(03114, 0);
	SaveEMEM(03115, 0);
	SaveEMEM(03116, 0);
	SaveEMEM(03117, 0);
	SaveEMEM(03120, 0);

	//AGSK
	DoubleToBuffer(BLOCKII.AGSK*3600.0*100.0, 28, iTemp, iTemp2);
	SaveEMEM(03371, iTemp);
	SaveEMEM(03372, iTemp2);

	//AZBIAS
	iTemp = SingleToBuffer(0.0 / 360.0, -1);
	SaveEMEM(03373, iTemp);

	//AOTAZ
	iTemp = SingleToBuffer(-60.0 / 360.0, -1, true);
	SaveEMEM(03404, iTemp);
	iTemp = SingleToBuffer(0.0 / 360.0, -1, true);
	SaveEMEM(03405, iTemp);
	iTemp = SingleToBuffer(60.0 / 360.0, -1, true);
	SaveEMEM(03406, iTemp);
	iTemp = SingleToBuffer(120.0 / 360.0, -1, true);
	SaveEMEM(03407, iTemp);
	iTemp = SingleToBuffer(-180.0 / 360.0, -1, true);
	SaveEMEM(03410, iTemp);
	iTemp = SingleToBuffer(-120.0 / 360.0, -1, true);
	SaveEMEM(03411, iTemp);
	//AOTEL
	iTemp = SingleToBuffer(45.0 / 360.0, -1);
	SaveEMEM(03412, iTemp);
	iTemp = SingleToBuffer(45.0 / 360.0, -1);
	SaveEMEM(03413, iTemp);
	iTemp = SingleToBuffer(45.0 / 360.0, -1);
	SaveEMEM(03414, iTemp);
	iTemp = SingleToBuffer(45.0 / 360.0, -1);
	SaveEMEM(03415, iTemp);
	iTemp = SingleToBuffer(45.0 / 360.0, -1);
	SaveEMEM(03416, iTemp);
	iTemp = SingleToBuffer(45.0 / 360.0, -1);
	SaveEMEM(03417, iTemp);

	//LRHMAX
	iTemp = SingleToBuffer(50000.0*0.3048, 14);
	SaveEMEM(03420, iTemp);
	//LRWH
	iTemp = SingleToBuffer(0.35, 0);
	SaveEMEM(03421, iTemp);
	//ZOOMTIME
	iTemp = SingleToBuffer(26.0*100.0, 14);
	SaveEMEM(03422, iTemp);
	//TENDBRAK
	iTemp = SingleToBuffer(62.0*100.0, 17);
	SaveEMEM(03423, iTemp);
	//TENDAPPR
	iTemp = SingleToBuffer(12.0*100.0, 17);
	SaveEMEM(03424, iTemp);
	//DELTTFAP
	iTemp = SingleToBuffer(BLOCKII.DELTTFAP*100.0, 17);
	SaveEMEM(03425, iTemp);
	//LEADTIME
	iTemp = SingleToBuffer(-2.2*100.0, 17);
	SaveEMEM(03426, iTemp);
	//LEADTIME
	iTemp = SingleToBuffer(62.0*100.0, 17);
	SaveEMEM(03427, iTemp);
	//RPCRTQSW
	iTemp = SingleToBuffer(-1.0, 1);
	SaveEMEM(03430, iTemp);
	//TNEWA
	SaveEMEM(03431, 020000);
	SaveEMEM(03432, 0);
}

void AGCPadloadGenerator::IMUCompensation(bool cmc, bool earlymodel)
{
	IMUBiasCompensationData *data;
	double pipapulse, MERUG;
	int pipascal, imuscal;

	const double RADTOPULSES = 2097152.0 / PI2;
	const double EARTHRATE = 15.0 / 3600.0*RAD;// 7.29209817e-5;
	const double GYROPULSE = RADTOPULSES* EARTHRATE*1e-5;// 0.00024272592; //gyro pulses/cs
	const double CMPIPAPULSE = 5.85; //cm/sec
	const double LMPIPAPULSE = 1.0; //cm/sec
	const double PPM = 1e-6; //parts per million
	const double GG = 9.7924; //m/sec^2

	//Variable scaling
	if (cmc)
	{
		data = &CMCDATA.IMUBiasCompensation;
		pipapulse = CMPIPAPULSE * 100.0;
		MERUG = GYROPULSE* CMPIPAPULSE / GG;
		imuscal = -3;

		if (earlymodel)
		{
			//Through Apollo 11
			pipascal = -8;
		}
		else
		{
			//Apollo 12 and later
			pipascal = -6;
		}
	}
	else
	{
		data = &LGCDATA.IMUBiasCompensation;
		pipapulse = LMPIPAPULSE * 100.0;
		MERUG = GYROPULSE* LMPIPAPULSE / GG*2.0; //TBD: Where does the 2.0 come from?
		imuscal = -5;

		if (earlymodel)
		{
			//Through Apollo 11
			pipascal = -5;
		}
		else
		{
			//Apollo 12 and later
			pipascal = -3;
		}
	}

	//PBIASX
	iTemp = SingleToBuffer(data->PBIASX / pipapulse, pipascal);
	SaveEMEM(01452, iTemp);
	//PIPASCFX
	iTemp = SingleToBuffer(data->PIPASCFX*PPM, -9);
	SaveEMEM(01453, iTemp);
	//PBIASY
	iTemp = SingleToBuffer(data->PBIASY / pipapulse, pipascal);
	SaveEMEM(01454, iTemp);
	//PIPASCFY
	iTemp = SingleToBuffer(data->PIPASCFY*PPM, -9);
	SaveEMEM(01455, iTemp);
	//PBIASZ
	iTemp = SingleToBuffer(data->PBIASZ / pipapulse, pipascal);
	SaveEMEM(01456, iTemp);
	//PIPASCFZ
	iTemp = SingleToBuffer(data->PIPASCFZ*PPM, -9);
	SaveEMEM(01457, iTemp);
	//NBDX
	iTemp = SingleToBuffer(data->NBDX*GYROPULSE, -5);
	SaveEMEM(01460, iTemp);
	//NBDY
	iTemp = SingleToBuffer(data->NBDY*GYROPULSE, -5);
	SaveEMEM(01461, iTemp);
	//NBDZ
	iTemp = SingleToBuffer(data->NBDZ*GYROPULSE, -5);
	SaveEMEM(01462, iTemp);
	//ADIAX
	iTemp = SingleToBuffer(data->ADIAX*MERUG, imuscal);
	SaveEMEM(01463, iTemp);
	//ADIAY
	iTemp = SingleToBuffer(data->ADIAY*MERUG, imuscal);
	SaveEMEM(01464, iTemp);
	//ADIAZ
	iTemp = SingleToBuffer(data->ADIAZ*MERUG, imuscal);
	SaveEMEM(01465, iTemp);
	//ADSRAX
	iTemp = SingleToBuffer(data->ADSRAX*MERUG, imuscal);
	SaveEMEM(01466, iTemp);
	//ADSRAY
	iTemp = SingleToBuffer(data->ADSRAY*MERUG, imuscal);
	SaveEMEM(01467, iTemp);
	//ADSRAZ
	iTemp = SingleToBuffer(data->ADSRAZ*MERUG, imuscal);
	SaveEMEM(01470, iTemp);
}

void AGCPadloadGenerator::R2Model(int address)
{
	if (BLOCKII.R2Model)
	{
		//E3J22R3M
		SaveEMEM(address, 012160);
		//E32C3IRM
		SaveEMEM(address + 1, 03363);
	}
	else
	{
		//E3J22R3M
		SaveEMEM(address, 0);
		//E32C3IRM
		SaveEMEM(address + 1, 0);
	}
}

void AGCPadloadGenerator::CMCDefaults(bool EarlyPIPABias, bool IsC108)
{
	//Contains padloads that never change their address in all of Colossus

	//IsC108: Comanche 108 specific address differences

	//CDUCHKWD (5 centiseconds)
	//Single precision erasable memory constant, program notation
	//"CDUCHKWD", scale factor B14, used to specify (if positive non -
	//zero) the number of centi-seconds delay before "MARKDIF" is performed
	//after receiving an optics mark button input. If the cell is zero
	//or negative, the delay is 0.01 second.
	iTemp = SingleToBuffer(CDUCHKWD * 100.0, 14);
	if (IsC108)
	{
		SaveEMEM(01346, iTemp);
	}
	else
	{
		SaveEMEM(01341, iTemp);
	}

	IMUCompensation(true, EarlyPIPABias);

	//WRENDPOS
	iTemp = SingleToBuffer(BLOCKII.WRENDPOS, 19);
	SaveEMEM(02000, iTemp);
	//WRENDVEL
	iTemp = SingleToBuffer(BLOCKII.WRENDVEL / 100.0, 0);
	SaveEMEM(02001, iTemp);

	//RMAX
	if (BLOCKII.RMAX > 0)
	{
		iTemp = SingleToBuffer(BLOCKII.RMAX*0.3048, 19);
		SaveEMEM(02002, iTemp);
	}
	else
	{
		SaveEMEM(02002, 077776);
	}

	//VMAX
	if (BLOCKII.VMAX > 0)
	{
		iTemp = SingleToBuffer(BLOCKII.VMAX*0.3048 / 100.0, 7);
		SaveEMEM(02003, iTemp);
	}
	else
	{
		SaveEMEM(02003, 077776);
	}

	//WORBPOS
	iTemp = SingleToBuffer(WORBPOS, 19);
	SaveEMEM(02004, iTemp);
	//WORBVEL
	iTemp = SingleToBuffer(WORBVEL / 100.0, 0);
	SaveEMEM(02005, iTemp);
	//S22WSUBL
	iTemp = SingleToBuffer(S22WSUBL, 19); //10km
	SaveEMEM(02006, iTemp);
	//RPVAR
	DoubleToBuffer(RPVAR, 28, iTemp, iTemp2); //4,000,000 m^2
	SaveEMEM(02007, iTemp);
	SaveEMEM(02010, iTemp2);

	//EMSALT
	DoubleToBuffer(EMSALT*0.3048, 29, iTemp, iTemp2);
	SaveEMEM(02017, iTemp);
	SaveEMEM(02020, iTemp2);

	//ATIGINC
	dTemp = 180.0*100.0;
	DoubleToBuffer(dTemp, 28, iTemp, iTemp2);
	SaveEMEM(02021, iTemp);
	SaveEMEM(02022, iTemp2);
	//PTIGINC
	DoubleToBuffer(dTemp, 28, iTemp, iTemp2);
	SaveEMEM(02023, iTemp);
	SaveEMEM(02024, iTemp2);

	//RLS

	RLS = r_from_latlong(LSLat*RAD, LSLng*RAD, LSAlt, BODY_MOON, 1);
	DoubleToBuffer(RLS.x, 27, iTemp, iTemp2);
	SaveEMEM(02025, iTemp);
	SaveEMEM(02026, iTemp2);
	DoubleToBuffer(RLS.y, 27, iTemp, iTemp2);
	SaveEMEM(02027, iTemp);
	SaveEMEM(02030, iTemp2);
	DoubleToBuffer(RLS.z, 27, iTemp, iTemp2);
	SaveEMEM(02031, iTemp);
	SaveEMEM(02032, iTemp2);

	//INTVAR
	dTemp = 196.0; //196 m^2
	iTemp = SingleToBuffer(dTemp, 15);
	SaveEMEM(02377, iTemp);

	//AZIMUTH
	dTemp = PadAzimuth / 360.0;
	DoubleToBuffer(dTemp, 0, iTemp, iTemp2);
	SaveEMEM(02400, iTemp);
	SaveEMEM(02401, iTemp2);

	//LATITUDE
	dTemp = PadLat / 360.0;
	DoubleToBuffer(dTemp, 0, iTemp, iTemp2);
	SaveEMEM(02402, iTemp);
	SaveEMEM(02403, iTemp2);

	//TAZEL
	dTemp = TAZEL[0] / 360.0;
	iTemp = SingleToBuffer(dTemp, -1);
	if (iTemp > 037777) iTemp++; //Two's complement
	SaveEMEM(02432, iTemp);
	dTemp = TAZEL[1] / 360.0;
	iTemp = SingleToBuffer(dTemp, -2);
	SaveEMEM(02433, iTemp);
	dTemp = TAZEL[2] / 360.0;
	iTemp = SingleToBuffer(dTemp, -1);
	if (iTemp > 037777) iTemp++; //Two's complement
	SaveEMEM(02434, iTemp);
	dTemp = TAZEL[3] / 360.0;
	iTemp = SingleToBuffer(dTemp, -2);
	SaveEMEM(02435, iTemp);

	//RVAR
	SaveEMEM(03002, 0);
	SaveEMEM(03003, 0);
	//RVARMIN
	TripleToBuffer(-BLOCKII.RVARMIN*pow(FT2M, 2), 40, iTemp, iTemp2, iTemp3);
	SaveEMEM(03004, iTemp);
	SaveEMEM(03005, iTemp2);
	SaveEMEM(03006, iTemp3);

	//LAT(SPL)
	dTemp = BLOCKII.LAT_SPL / 360.0;
	DoubleToBuffer(dTemp, 0, iTemp, iTemp2);
	SaveEMEM(03400, iTemp);
	SaveEMEM(03401, iTemp2);
	//LNG(SPL)
	dTemp = BLOCKII.LNG_SPL / 360.0;
	DoubleToBuffer(dTemp, 0, iTemp, iTemp2);
	SaveEMEM(03402, iTemp);
	SaveEMEM(03403, iTemp2);
}

void AGCPadloadGenerator::Colossus237_249_Defaults(bool Is249)
{
	CMCDefaults(true);

	if (Is249)
	{
		//EMDOT
		dTemp = SPS_THRUST / SPS_ISP;
		iTemp = SingleToBuffer(dTemp / 100.0, 3);
		SaveEMEM(0110, iTemp);
	}

	//PIPTIME
	double PIPTIME = (LaunchMJD - PrelaunchMJD)*8.64e6;
	DoubleToBuffer(PIPTIME, 28, iTemp, iTemp2);
	SaveEMEM(01204, iTemp);
	SaveEMEM(01205, iTemp2);

	//PADLONG
	dTemp = PadLong / 360.0;
	DoubleToBuffer(dTemp, 0, iTemp, iTemp2);
	SaveEMEM(01263, iTemp);
	SaveEMEM(01264, iTemp2);

	//PGNCSALT
	dTemp = PadAlt;
	DoubleToBuffer(dTemp, 29, iTemp, iTemp2);
	SaveEMEM(01272, iTemp);
	SaveEMEM(01273, iTemp2);

	//RTED1
	DoubleToBuffer(RTED1, 3, iTemp, iTemp2);
	SaveEMEM(01351, iTemp);
	SaveEMEM(01352, iTemp2);
	//DVTHRESH
	iTemp = SingleToBuffer(DVTHRESH / 100.0, -2);
	SaveEMEM(01353, iTemp);
	//HORIZALT
	DoubleToBuffer(HORIZALT, 29, iTemp, iTemp2);
	SaveEMEM(01354, iTemp);
	SaveEMEM(01355, iTemp2);
	//ALTVAR
	iTemp = SingleToBuffer(ALTVAR, -16);
	SaveEMEM(01356, iTemp);

	//LAUNCHAZ
	dTemp = LaunchAzimuth / 360.0;
	DoubleToBuffer(dTemp, 0, iTemp, iTemp2);
	SaveEMEM(02633, iTemp);
	SaveEMEM(02634, iTemp2);

	//WMIDPOS
	dTemp = CMCDATA.WMIDPOS*FT2M / sqrt(3.0);
	iTemp = SingleToBuffer(dTemp, 19);
	SaveEMEM(03000, iTemp);
	//WMIDVEL
	dTemp = CMCDATA.WMIDVEL*FT2M / 100.0 / sqrt(3.0);
	iTemp = SingleToBuffer(dTemp, 0);
	SaveEMEM(03001, iTemp);

	//LADPAD
	dTemp = 0.3;
	iTemp = SingleToBuffer(dTemp, 0);
	SaveEMEM(03007, iTemp);
	//LODPAD
	dTemp = 0.18;
	iTemp = SingleToBuffer(dTemp, 0);
	SaveEMEM(03010, iTemp);
	//ALFAPAD
	dTemp = -18.51 / 360.0;
	iTemp = SingleToBuffer(dTemp, -1);
	SaveEMEM(03011, iTemp);

	//ESTROKER
	SaveEMEM(03012, 05);
	//EKPRIME
	SaveEMEM(03013, 0123);
	SaveEMEM(03014, 077);
	//ETDECAY
	dTemp = 0.6*100.0;
	iTemp = SingleToBuffer(dTemp, 14);
	SaveEMEM(03015, iTemp);
	//EKTLX
	SaveEMEM(03016, 017433);
	SaveEMEM(03017, 04222);
	//ETVCDT/2
	SaveEMEM(03020, 02);
	//ETSWITCH
	SaveEMEM(03021, 06);
	//ECORFRAC
	SaveEMEM(03022, 010000);
	//EREPFRAC
	SaveEMEM(03023, 01000);
	SaveEMEM(03024, 0232);

	//PACTOFF
	SaveEMEM(03025, 077676);
	//YACTOFF
	SaveEMEM(03026, 070);

	//AP0
	SaveEMEM(03027, 010000);
	//AP1
	SaveEMEM(03030, 075024);
	SaveEMEM(03031, 047207);
	//AP2
	SaveEMEM(03032, 01066);
	SaveEMEM(03033, 033013);
	//AP3
	SaveEMEM(03034, 077705);
	SaveEMEM(03035, 061340);
	//BP1
	SaveEMEM(03036, 057377);
	SaveEMEM(03037, 074664);
	//BP2
	SaveEMEM(03040, 014414);
	SaveEMEM(03041, 036243);
	//BP3
	SaveEMEM(03042, 074624);
	SaveEMEM(03043, 073117);

	//DAPDATR1
	SaveEMEM(03066, 031102);
	//DAPDATR2
	SaveEMEM(03067, 01111);

	//LEMMASS
	iTemp = SingleToBuffer(BLOCKII.LMMass*LBS2KG, 16);
	SaveEMEM(03073, iTemp);
	//CSMMASS
	iTemp = SingleToBuffer(BLOCKII.CSMMass*LBS2KG, 16);
	SaveEMEM(03074, iTemp);

	//POLYNUM
	SaveEMEM(03261, 05);
	SaveEMEM(03262, 077777);
	SaveEMEM(03263, 047314);
	SaveEMEM(03264, 026);
	SaveEMEM(03265, 021213);
	SaveEMEM(03266, 0573);
	SaveEMEM(03267, 011306);
	SaveEMEM(03270, 077256);
	SaveEMEM(03271, 063275);
	SaveEMEM(03272, 0);
	SaveEMEM(03273, 0);
	SaveEMEM(03274, 0);
	SaveEMEM(03275, 0);
	SaveEMEM(03276, 0);
	SaveEMEM(03277, 0);

	//SATRLRT
	SaveEMEM(03300, 0);
	SaveEMEM(03301, 012144);
	//RPSTART
	dTemp = 11.85*100.0;
	iTemp = SingleToBuffer(dTemp, 14);
	SaveEMEM(03302, iTemp);
	//POLYSTOP
	dTemp = -147.25*100.0;
	iTemp = SingleToBuffer(dTemp, 14);
	SaveEMEM(03303, iTemp);

	//ECSTEER
	SaveEMEM(03424, 010000);
}

void AGCPadloadGenerator::Comanche45Padload(bool IsR2)
{
	CMCDefaults(true);

	//FLAGWRD1
	SaveEMEM(075, 0);
	//FLAGWRD3
	SaveEMEM(077, 0);
	//FLAGWRD8
	SaveEMEM(0104, 0);
	//EMDOT
	dTemp = CMCDATA.EMDOT*LBS2KG;
	iTemp = SingleToBuffer(dTemp / 100.0, 3);
	SaveEMEM(0110, iTemp);
	//DUMPCNT
	SaveEMEM(0333, 010000);
	//PIPTIME
	double PIPTIME = (LaunchMJD - PrelaunchMJD)*8.64e6;
	DoubleToBuffer(PIPTIME, 28, iTemp, iTemp2);
	SaveEMEM(01204, iTemp);
	SaveEMEM(01205, iTemp2);
	//PADLONG
	dTemp = PadLong / 360.0;
	DoubleToBuffer(dTemp, 0, iTemp, iTemp2);
	SaveEMEM(01263, iTemp);
	SaveEMEM(01264, iTemp2);
	//PGNCSALT
	dTemp = PadAlt;
	DoubleToBuffer(dTemp, 29, iTemp, iTemp2);
	SaveEMEM(01272, iTemp);
	SaveEMEM(01273, iTemp2);

	//RTED1
	DoubleToBuffer(RTED1, 3, iTemp, iTemp2);
	SaveEMEM(01351, iTemp);
	SaveEMEM(01352, iTemp2);

	//DVTHRESH
	iTemp = SingleToBuffer(DVTHRESH / 100.0, -2);
	SaveEMEM(01353, iTemp);

	//HORIZALT
	DoubleToBuffer(HORIZALT, 29, iTemp, iTemp2);
	SaveEMEM(01354, iTemp);
	SaveEMEM(01355, iTemp2);
	//ALTVAR
	iTemp = SingleToBuffer(ALTVAR, -16);
	SaveEMEM(01356, iTemp);

	//GCOMPSW
	SaveEMEM(01477, 0);

	//EK1VAL
	DoubleToBuffer(CMCDATA.EK1VAL*LBF2N / 100.0, 23, iTemp, iTemp2);
	SaveEMEM(01767, iTemp);
	SaveEMEM(01770, iTemp2);
	//FANG
	DoubleToBuffer(CMCDATA.FANG*LBF2N / pow(100.0, 2), 7, iTemp, iTemp2);
	SaveEMEM(01771, iTemp);
	SaveEMEM(01772, iTemp2);

	if (IsR2)
	{
		R2Model(01773);
	}

	//LAUNCHAZ
	dTemp = LaunchAzimuth / 360.0;
	DoubleToBuffer(dTemp, 0, iTemp, iTemp2);
	SaveEMEM(02633, iTemp);
	SaveEMEM(02634, iTemp2);

	//WMIDPOS
	dTemp = CMCDATA.WMIDPOS*FT2M / sqrt(3.0);
	iTemp = SingleToBuffer(dTemp, 19);
	SaveEMEM(03000, iTemp);
	//WMIDVEL
	dTemp = CMCDATA.WMIDVEL*FT2M / 100.0 / sqrt(3.0);
	iTemp = SingleToBuffer(dTemp, 0);
	SaveEMEM(03001, iTemp);

	//LADPAD
	dTemp = BLOCKII.LADPAD;
	iTemp = SingleToBuffer(dTemp, 0);
	SaveEMEM(03007, iTemp);
	//LODPAD
	dTemp = BLOCKII.LODPAD;
	iTemp = SingleToBuffer(dTemp, 0);
	SaveEMEM(03010, iTemp);
	//ALFAPAD
	dTemp = BLOCKII.ALFAPAD / 360.0;
	iTemp = SingleToBuffer(dTemp, -1);
	SaveEMEM(03011, iTemp);

	//ETDECAY
	dTemp = 0.6*100.0;
	iTemp = SingleToBuffer(dTemp, 14);
	SaveEMEM(03012, iTemp);
	//ESTROKER
	SaveEMEM(03013, 05);
	//EKPRIME
	SaveEMEM(03014, 0123);
	SaveEMEM(03015, 0175);

	//EKTLX
	SaveEMEM(03016, 037065);
	SaveEMEM(03017, 02245);
	SaveEMEM(03020, 0156);
	//EREPFRAC
	SaveEMEM(03021, 01000);
	SaveEMEM(03022, 0232);
	//PACTOFF
	iTemp = SingleToBuffer(BLOCKII.PACTOFF / 0.023725, 14);
	SaveEMEM(03023, iTemp);
	//YACTOFF
	iTemp = SingleToBuffer(BLOCKII.YACTOFF / 0.023725, 14);
	SaveEMEM(03024, iTemp);
	//HBN10
	SaveEMEM(03025, 037777);
	//HBN11/2
	SaveEMEM(03026, 0);
	//HBN12
	SaveEMEM(03027, 0);
	//HBD11/2
	SaveEMEM(03030, 054360);
	//HBD12
	SaveEMEM(03031, 021075);
	//HBN20
	SaveEMEM(03032, 037777);
	//HBN21/2
	SaveEMEM(03033, 060465);
	//HBN22
	SaveEMEM(03034, 0);
	//HBD21/2
	SaveEMEM(03035, 054360);
	//HBD22
	SaveEMEM(03036, 021075);
	//HBN30
	SaveEMEM(03037, 037777);
	//HBN31/2
	SaveEMEM(03040, 057142);
	//HBN32
	SaveEMEM(03041, 033106);
	//HBD31/2
	SaveEMEM(03042, 050741);
	//HBD32
	SaveEMEM(03043, 031162);
	//DAPDATR1
	SaveEMEM(03066, 031102);
	//DAPDATR2
	SaveEMEM(03067, 01111);
	//LEMMASS
	iTemp = SingleToBuffer(BLOCKII.LMMass*LBS2KG, 16);
	SaveEMEM(03073, iTemp);
	//CSMMASS
	iTemp = SingleToBuffer(BLOCKII.CSMMass*LBS2KG, 16);
	SaveEMEM(03074, iTemp);
	//POLYNUM
	SaveEMEM(03261, 05);
	SaveEMEM(03262, 0);
	SaveEMEM(03263, 014107);
	SaveEMEM(03264, 010);
	SaveEMEM(03265, 024775);
	SaveEMEM(03266, 0536);
	SaveEMEM(03267, 010550);
	SaveEMEM(03270, 0746);
	SaveEMEM(03271, 027545);
	SaveEMEM(03272, 072636);
	SaveEMEM(03273, 060425);
	SaveEMEM(03274, 05754);
	SaveEMEM(03275, 036030);
	SaveEMEM(03276, 075607);
	SaveEMEM(03277, 053465);
	//SATRLRT
	SaveEMEM(03300, 0);
	SaveEMEM(03301, 016441);
	//RPSTART
	dTemp = 11.85*100.0;
	iTemp = SingleToBuffer(dTemp, 14);
	SaveEMEM(03302, iTemp);
	//POLYSTOP
	dTemp = -147.25*100.0;
	iTemp = SingleToBuffer(dTemp, 14);
	SaveEMEM(03303, iTemp);
	//SATRATE
	SaveEMEM(03321, 0);
	SaveEMEM(03322, 0344);
	SaveEMEM(03323, 077433);
	SaveEMEM(03324, 0);
	//SATSCALE
	SaveEMEM(03331, 010000); //0.3 V/DEG
	//P37RANGE
	dTemp = BLOCKII.P37RANGE / 21600.0;
	iTemp = SingleToBuffer(dTemp, 0);
	SaveEMEM(03376, iTemp);

	//ECSTEER
	SaveEMEM(03424, 010000);
}

void AGCPadloadGenerator::Comanche55Padload()
{
	CMCDefaults(true);

	//FLAGWRD1
	SaveEMEM(075, 0);
	//FLAGWRD3
	SaveEMEM(077, 0);
	//FLAGWRD8
	SaveEMEM(0104, 0);
	//EMDOT
	dTemp = CMCDATA.EMDOT*LBS2KG;
	iTemp = SingleToBuffer(dTemp / 100.0, 3);
	SaveEMEM(0110, iTemp);
	//DUMPCNT
	SaveEMEM(0333, 010000);
	//PIPTIME
	double PIPTIME = (LaunchMJD - PrelaunchMJD)*8.64e6;
	DoubleToBuffer(PIPTIME, 28, iTemp, iTemp2);
	SaveEMEM(01204, iTemp);
	SaveEMEM(01205, iTemp2);
	//PADLONG
	dTemp = PadLong / 360.0;
	DoubleToBuffer(dTemp, 0, iTemp, iTemp2);
	SaveEMEM(01263, iTemp);
	SaveEMEM(01264, iTemp2);
	//PGNCSALT
	dTemp = PadAlt;
	DoubleToBuffer(dTemp, 29, iTemp, iTemp2);
	SaveEMEM(01272, iTemp);
	SaveEMEM(01273, iTemp2);

	//RTED1
	DoubleToBuffer(RTED1, 3, iTemp, iTemp2);
	SaveEMEM(01351, iTemp);
	SaveEMEM(01352, iTemp2);

	//DVTHRESH
	iTemp = SingleToBuffer(DVTHRESH / 100.0, -2);
	SaveEMEM(01353, iTemp);

	//HORIZALT
	DoubleToBuffer(HORIZALT, 29, iTemp, iTemp2);
	SaveEMEM(01354, iTemp);
	SaveEMEM(01355, iTemp2);
	//ALTVAR
	iTemp = SingleToBuffer(ALTVAR, -16);
	SaveEMEM(01356, iTemp);

	//GCOMPSW
	SaveEMEM(01477, 0);

	//EK1VAL
	DoubleToBuffer(CMCDATA.EK1VAL*LBF2N / 100.0, 23, iTemp, iTemp2);
	SaveEMEM(01767, iTemp);
	SaveEMEM(01770, iTemp2);
	//EK2VAL
	DoubleToBuffer(CMCDATA.EK2VAL*LBF2N / 100.0, 23, iTemp, iTemp2);
	SaveEMEM(01771, iTemp);
	SaveEMEM(01772, iTemp2);
	//EK3VAL
	iTemp = SingleToBuffer(CMCDATA.EK3VAL*LBF2N / pow(100.0, 2), 9);
	SaveEMEM(01773, iTemp);
	//FANG
	DoubleToBuffer(CMCDATA.FANG*LBF2N / pow(100.0, 2), 7, iTemp, iTemp2);
	SaveEMEM(01774, iTemp);
	SaveEMEM(01775, iTemp2);

	R2Model(01776);

	//LAUNCHAZ
	dTemp = LaunchAzimuth / 360.0;
	DoubleToBuffer(dTemp, 0, iTemp, iTemp2);
	SaveEMEM(02633, iTemp);
	SaveEMEM(02634, iTemp2);

	//WMIDPOS
	dTemp = CMCDATA.WMIDPOS*FT2M;
	iTemp = SingleToBuffer(dTemp, 19);
	SaveEMEM(03000, iTemp);
	//WMIDVEL
	dTemp = CMCDATA.WMIDVEL*FT2M / 100.0;
	iTemp = SingleToBuffer(dTemp, 0);
	SaveEMEM(03001, iTemp);

	//LADPAD
	dTemp = BLOCKII.LADPAD;
	iTemp = SingleToBuffer(dTemp, 0);
	SaveEMEM(03007, iTemp);
	//LODPAD
	dTemp = BLOCKII.LODPAD;
	iTemp = SingleToBuffer(dTemp, 0);
	SaveEMEM(03010, iTemp);
	//ALFAPAD
	dTemp = BLOCKII.ALFAPAD / 360.0;
	iTemp = SingleToBuffer(dTemp, -1);
	SaveEMEM(03011, iTemp);

	//ETDECAY
	dTemp = 0.6*100.0;
	iTemp = SingleToBuffer(dTemp, 14);
	SaveEMEM(03012, iTemp);
	//ESTROKER
	SaveEMEM(03013, 05);
	//EKPRIME
	SaveEMEM(03014, 0123);
	SaveEMEM(03015, 0175);

	//EKTLX
	SaveEMEM(03016, 037065);
	SaveEMEM(03017, 02245);
	SaveEMEM(03020, 0156);
	//EREPFRAC
	SaveEMEM(03021, 01000);
	SaveEMEM(03022, 0232);
	//PACTOFF
	iTemp = SingleToBuffer(BLOCKII.PACTOFF / 0.023725, 14);
	SaveEMEM(03023, iTemp);
	//YACTOFF
	iTemp = SingleToBuffer(BLOCKII.YACTOFF / 0.023725, 14);
	SaveEMEM(03024, iTemp);
	//HBN10
	SaveEMEM(03025, 037777);
	//HBN11/2
	SaveEMEM(03026, 0);
	//HBN12
	SaveEMEM(03027, 0);
	//HBD11/2
	SaveEMEM(03030, 054360);
	//HBD12
	SaveEMEM(03031, 021075);
	//HBN20
	SaveEMEM(03032, 037777);
	//HBN21/2
	SaveEMEM(03033, 060465);
	//HBN22
	SaveEMEM(03034, 0);
	//HBD21/2
	SaveEMEM(03035, 054360);
	//HBD22
	SaveEMEM(03036, 021075);
	//HBN30
	SaveEMEM(03037, 037777);
	//HBN31/2
	SaveEMEM(03040, 057142);
	//HBN32
	SaveEMEM(03041, 033106);
	//HBD31/2
	SaveEMEM(03042, 050741);
	//HBD32
	SaveEMEM(03043, 031162);
	//DAPDATR1
	SaveEMEM(03066, 031102);
	//DAPDATR2
	SaveEMEM(03067, 01111);
	//LEMMASS
	iTemp = SingleToBuffer(BLOCKII.LMMass*LBS2KG, 16);
	SaveEMEM(03073, iTemp);
	//CSMMASS
	iTemp = SingleToBuffer(BLOCKII.CSMMass*LBS2KG, 16);
	SaveEMEM(03074, iTemp);
	//POLYNUM
	SavePOLYNUM(03261);
	//SATRLRT
	SaveEMEM(03300, 0);
	SaveEMEM(03301, 016441);
	//RPSTART
	dTemp = BLOCKII.RPSTART*100.0;
	iTemp = SingleToBuffer(dTemp, 14);
	SaveEMEM(03302, iTemp);
	//POLYSTOP
	dTemp = -BLOCKII.POLYSTOP*100.0;
	iTemp = SingleToBuffer(dTemp, 14);
	SaveEMEM(03303, iTemp);
	//SATRATE
	SaveEMEM(03321, 0);
	SaveEMEM(03322, 0344);
	SaveEMEM(03323, 077433);
	SaveEMEM(03324, 0);
	//SATSCALE
	SaveEMEM(03331, 010000); //0.3 V/DEG
	//P37RANGE
	dTemp = BLOCKII.P37RANGE / 21600.0;
	iTemp = SingleToBuffer(dTemp, 0);
	SaveEMEM(03376, iTemp);

	//ECSTEER
	SaveEMEM(03424, 010000);
}

void AGCPadloadGenerator::Comanche67Padload()
{
	CMCDefaults(false);

	//FLAGWRD1
	SaveEMEM(075, 0);
	//FLAGWRD3
	SaveEMEM(077, 0);
	//FLAGWRD8
	//Bit 8 set: SURFFLAG. LM on lunar surface
	SaveEMEM(0104, 0200);
	//EMDOT
	dTemp = CMCDATA.EMDOT*LBS2KG;
	iTemp = SingleToBuffer(dTemp / 100.0, 3);
	SaveEMEM(0110, iTemp);
	//PIPTIME
	double PIPTIME = (LaunchMJD - PrelaunchMJD)*8.64e6;
	DoubleToBuffer(PIPTIME, 28, iTemp, iTemp2);
	SaveEMEM(01204, iTemp);
	SaveEMEM(01205, iTemp2);
	//PADLONG
	dTemp = PadLong / 360.0;
	DoubleToBuffer(dTemp, 0, iTemp, iTemp2);
	SaveEMEM(01263, iTemp);
	SaveEMEM(01264, iTemp2);
	//PGNCSALT
	dTemp = PadAlt;
	DoubleToBuffer(dTemp, 29, iTemp, iTemp2);
	SaveEMEM(01272, iTemp);
	SaveEMEM(01273, iTemp2);

	//RTED1
	DoubleToBuffer(RTED1, 3, iTemp, iTemp2);
	SaveEMEM(01345, iTemp);
	SaveEMEM(01346, iTemp2);

	//DVTHRESH
	iTemp = SingleToBuffer(DVTHRESH / 100.0, -2);
	SaveEMEM(01347, iTemp);

	//HORIZALT
	DoubleToBuffer(HORIZALT, 29, iTemp, iTemp2);
	SaveEMEM(01350, iTemp);
	SaveEMEM(01351, iTemp2);
	//ALTVAR
	iTemp = SingleToBuffer(ALTVAR, -16);
	SaveEMEM(01352, iTemp);

	//GCOMPSW
	SaveEMEM(01477, 0);

	//EK1VAL
	DoubleToBuffer(CMCDATA.EK1VAL*LBF2N / 100.0, 23, iTemp, iTemp2);
	SaveEMEM(01767, iTemp);
	SaveEMEM(01770, iTemp2);
	//EK2VAL
	DoubleToBuffer(CMCDATA.EK2VAL*LBF2N / 100.0, 23, iTemp, iTemp2);
	SaveEMEM(01771, iTemp);
	SaveEMEM(01772, iTemp2);
	//EK3VAL
	iTemp = SingleToBuffer(CMCDATA.EK3VAL*LBF2N / pow(100.0, 2), 9);
	SaveEMEM(01773, iTemp);
	//FANG
	DoubleToBuffer(CMCDATA.FANG*LBF2N / pow(100.0, 2), 7, iTemp, iTemp2);
	SaveEMEM(01774, iTemp);
	SaveEMEM(01775, iTemp2);

	R2Model(01776);

	//LAUNCHAZ
	dTemp = LaunchAzimuth / 360.0;
	DoubleToBuffer(dTemp, 0, iTemp, iTemp2);
	SaveEMEM(02633, iTemp);
	SaveEMEM(02634, iTemp2);

	//WMIDPOS
	dTemp = CMCDATA.WMIDPOS*FT2M;
	iTemp = SingleToBuffer(dTemp, 19);
	SaveEMEM(03000, iTemp);
	//WMIDVEL
	dTemp = CMCDATA.WMIDVEL*FT2M / 100.0;
	iTemp = SingleToBuffer(dTemp, 0);
	SaveEMEM(03001, iTemp);

	//LADPAD
	dTemp = 0.27;
	iTemp = SingleToBuffer(dTemp, 0);
	SaveEMEM(03007, iTemp);
	//LODPAD
	dTemp = 0.207;
	iTemp = SingleToBuffer(dTemp, 0);
	SaveEMEM(03010, iTemp);
	//ALFAPAD
	dTemp = -20.5 / 360.0;
	iTemp = SingleToBuffer(dTemp, -1);
	SaveEMEM(03011, iTemp);

	//ETDECAY
	dTemp = 0.6*100.0;
	iTemp = SingleToBuffer(dTemp, 14);
	SaveEMEM(03012, iTemp);
	//ESTROKER
	SaveEMEM(03013, 05);
	//EKPRIME
	SaveEMEM(03014, 0123);
	SaveEMEM(03015, 0175);

	//EKTLX
	SaveEMEM(03016, 017433);
	SaveEMEM(03017, 04500);
	SaveEMEM(03020, 0334);
	//EREPFRAC
	SaveEMEM(03021, 01000);
	SaveEMEM(03022, 0232);
	//PACTOFF
	iTemp = SingleToBuffer(BLOCKII.PACTOFF / 0.023725, 14);
	SaveEMEM(03023, iTemp);
	//YACTOFF
	iTemp = SingleToBuffer(BLOCKII.YACTOFF / 0.023725, 14);
	SaveEMEM(03024, iTemp);
	//HBN10
	SaveEMEM(03025, 037777);
	//HBN11/2
	SaveEMEM(03026, 0);
	//HBN12
	SaveEMEM(03027, 0);
	//HBD11/2
	SaveEMEM(03030, 054360);
	//HBD12
	SaveEMEM(03031, 021075);
	//HBN20
	SaveEMEM(03032, 037777);
	//HBN21/2
	SaveEMEM(03033, 060465);
	//HBN22
	SaveEMEM(03034, 0);
	//HBD21/2
	SaveEMEM(03035, 054360);
	//HBD22
	SaveEMEM(03036, 021075);
	//HBN30
	SaveEMEM(03037, 037777);
	//HBN31/2
	SaveEMEM(03040, 057142);
	//HBN32
	SaveEMEM(03041, 033106);
	//HBD31/2
	SaveEMEM(03042, 050741);
	//HBD32
	SaveEMEM(03043, 031162);
	//DAPDATR1
	SaveEMEM(03066, 031102);
	//DAPDATR2
	SaveEMEM(03067, 01111);
	//LEMMASS
	iTemp = SingleToBuffer(BLOCKII.LMMass*LBS2KG, 16);
	SaveEMEM(03073, iTemp);
	//CSMMASS
	iTemp = SingleToBuffer(BLOCKII.CSMMass*LBS2KG, 16);
	SaveEMEM(03074, iTemp);
	//POLYNUM
	SavePOLYNUM(03261);
	//SATRLRT
	SaveEMEM(03300, 0);
	SaveEMEM(03301, 016441);
	//RPSTART
	dTemp = BLOCKII.RPSTART*100.0;
	iTemp = SingleToBuffer(dTemp, 14);
	SaveEMEM(03302, iTemp);
	//POLYSTOP
	dTemp = -BLOCKII.POLYSTOP*100.0;
	iTemp = SingleToBuffer(dTemp, 14);
	SaveEMEM(03303, iTemp);
	//SATRATE
	SaveEMEM(03321, 0);
	SaveEMEM(03322, 0344);
	SaveEMEM(03323, 077433);
	SaveEMEM(03324, 0);
	//SATSCALE
	SaveEMEM(03331, 010000); //0.3 V/DEG
	//P37RANGE
	dTemp = BLOCKII.P37RANGE / 21600.0;
	iTemp = SingleToBuffer(dTemp, 0);
	SaveEMEM(03376, iTemp);

	//ECSTEER
	SaveEMEM(03424, 010000);

	//TB6JOB
	SaveEMEM(03513, 04);
	SaveEMEM(03514, 06);
	SaveEMEM(03515, 031413);
	SaveEMEM(03516, 053140);
	SaveEMEM(03517, 06);
	SaveEMEM(03520, 040025);
	SaveEMEM(03521, 021140);
	SaveEMEM(03522, 05357);
	SaveEMEM(03523, 01550);
	SaveEMEM(03524, 010067);
	SaveEMEM(03525, 011573);
	SaveEMEM(03526, 04106);

	SaveEMEM(03527, 035055);
	SaveEMEM(03530, 05251);
	SaveEMEM(03531, 01543);
	SaveEMEM(03532, 06);
	SaveEMEM(03533, 041413);
	SaveEMEM(03534, 053046);
	SaveEMEM(03535, 01574);
	SaveEMEM(03537, 021046);
	SaveEMEM(03540, 033300);
	SaveEMEM(03541, 04676);
	SaveEMEM(03542, 020707);
	SaveEMEM(03543, 037663);
	SaveEMEM(03544, 05150);

	SaveEMEM(03545, 01525);
	SaveEMEM(03546, 010067);
	SaveEMEM(03547, 05340);
	SaveEMEM(03550, 035017);
	SaveEMEM(03551, 06);
	SaveEMEM(03552, 05012);
	SaveEMEM(03553, 035031);
	SaveEMEM(03554, 06);
	SaveEMEM(03555, 05011);
	SaveEMEM(03556, 06);
	SaveEMEM(03557, 030025);
	SaveEMEM(03560, 021337);
	SaveEMEM(03561, 05303);

	SaveEMEM(03562, 01750);
	SaveEMEM(03563, 045017);
	SaveEMEM(03564, 06);
	SaveEMEM(03565, 03012);
	SaveEMEM(03566, 045031);
	SaveEMEM(03567, 06);
	SaveEMEM(03570, 03011);
	SaveEMEM(03571, 025573);
	SaveEMEM(03572, 05340);
	SaveEMEM(03573, 0);
	SaveEMEM(03574, 06);
	SaveEMEM(03575, 030025);
	SaveEMEM(03576, 01537);
}

void AGCPadloadGenerator::Comanche72Padload()
{
	CMCDefaults(false);

	//FLAGWRD1
	SaveEMEM(075, 0);
	//FLAGWRD3
	SaveEMEM(077, 0);
	//FLAGWRD8
	//Bit 8 set: SURFFLAG. LM on lunar surface
	SaveEMEM(0104, 0200);
	//EMDOT
	dTemp = CMCDATA.EMDOT*LBS2KG;
	iTemp = SingleToBuffer(dTemp / 100.0, 3);
	SaveEMEM(0110, iTemp);
	//PIPTIME
	double PIPTIME = (LaunchMJD - PrelaunchMJD)*8.64e6;
	DoubleToBuffer(PIPTIME, 28, iTemp, iTemp2);
	SaveEMEM(01204, iTemp);
	SaveEMEM(01205, iTemp2);
	//PADLONG
	dTemp = PadLong / 360.0;
	DoubleToBuffer(dTemp, 0, iTemp, iTemp2);
	SaveEMEM(01263, iTemp);
	SaveEMEM(01264, iTemp2);
	//PGNCSALT
	dTemp = PadAlt;
	DoubleToBuffer(dTemp, 29, iTemp, iTemp2);
	SaveEMEM(01272, iTemp);
	SaveEMEM(01273, iTemp2);
	//RTED1. First coefficient defining high speed V-gamma target line polynomial
	DoubleToBuffer(RTED1, 3, iTemp, iTemp2);
	SaveEMEM(01345, iTemp);
	SaveEMEM(01346, iTemp2);
	//DVTHRESH
	iTemp = SingleToBuffer(DVTHRESH / 100.0, -2);
	SaveEMEM(01347, iTemp);
	//HORIZALT
	DoubleToBuffer(HORIZALT, 29, iTemp, iTemp2);
	SaveEMEM(01350, iTemp);
	SaveEMEM(01351, iTemp2);
	//ALTVAR
	iTemp = SingleToBuffer(ALTVAR, -16);
	SaveEMEM(01352, iTemp);
	//EK1VAL
	DoubleToBuffer(CMCDATA.EK1VAL*LBF2N / 100.0, 23, iTemp, iTemp2);
	SaveEMEM(01767, iTemp);
	SaveEMEM(01770, iTemp2);
	//EK2VAL
	DoubleToBuffer(CMCDATA.EK2VAL*LBF2N / 100.0, 23, iTemp, iTemp2);
	SaveEMEM(01771, iTemp);
	SaveEMEM(01772, iTemp2);
	//EK3VAL
	iTemp = SingleToBuffer(CMCDATA.EK3VAL*LBF2N / pow(100.0, 2), 9);
	SaveEMEM(01773, iTemp);
	//FANG
	DoubleToBuffer(CMCDATA.FANG*LBF2N / pow(100.0, 2), 7, iTemp, iTemp2);
	SaveEMEM(01774, iTemp);
	SaveEMEM(01775, iTemp2);

	R2Model(01776);

	//LAUNCHAZ
	dTemp = LaunchAzimuth / 360.0;
	DoubleToBuffer(dTemp, 0, iTemp, iTemp2);
	SaveEMEM(02633, iTemp);
	SaveEMEM(02634, iTemp2);

	//WMIDPOS
	dTemp = CMCDATA.WMIDPOS*FT2M;
	iTemp = SingleToBuffer(dTemp, 19);
	SaveEMEM(03000, iTemp);
	//WMIDVEL
	dTemp = CMCDATA.WMIDVEL*FT2M / 100.0;
	iTemp = SingleToBuffer(dTemp, 0);
	SaveEMEM(03001, iTemp);

	//LADPAD
	dTemp = 0.3;
	iTemp = SingleToBuffer(dTemp, 0);
	SaveEMEM(03007, iTemp);
	//LODPAD
	dTemp = 0.18;
	iTemp = SingleToBuffer(dTemp, 0);
	SaveEMEM(03010, iTemp);
	//ALFAPAD
	dTemp = -21.49 / 360.0;
	iTemp = SingleToBuffer(dTemp, -1);
	SaveEMEM(03011, iTemp);
	//ETDECAY
	dTemp = 0.6*100.0;
	iTemp = SingleToBuffer(dTemp, 14);
	SaveEMEM(03012, iTemp);
	//ESTROKER
	SaveEMEM(03013, 05);
	//EKPRIME
	SaveEMEM(03014, 0123);
	SaveEMEM(03015, 0175);
	//EKTLX
	SaveEMEM(03016, 017433);
	SaveEMEM(03017, 04500);
	SaveEMEM(03020, 0334);
	//EREPFRAC
	SaveEMEM(03021, 01000);
	SaveEMEM(03022, 0232);
	//PACTOFF
	iTemp = SingleToBuffer(BLOCKII.PACTOFF / 0.023725, 14);
	SaveEMEM(03023, iTemp);
	//YACTOFF
	iTemp = SingleToBuffer(BLOCKII.YACTOFF / 0.023725, 14);
	SaveEMEM(03024, iTemp);
	//HBN10
	SaveEMEM(03025, 037777);
	//HBN11/2
	SaveEMEM(03026, 0);
	//HBN12
	SaveEMEM(03027, 0);
	//HBD11/2
	SaveEMEM(03030, 054360);
	//HBD12
	SaveEMEM(03031, 021075);
	//HBN20
	SaveEMEM(03032, 037777);
	//HBN21/2
	SaveEMEM(03033, 060465);
	//HBN22
	SaveEMEM(03034, 0);
	//HBD21/2
	SaveEMEM(03035, 054360);
	//HBD22
	SaveEMEM(03036, 021075);
	//HBN30
	SaveEMEM(03037, 037777);
	//HBN31/2
	SaveEMEM(03040, 057142);
	//HBN32
	SaveEMEM(03041, 033106);
	//HBD31/2
	SaveEMEM(03042, 050741);
	//HBD32
	SaveEMEM(03043, 031162);
	//DAPDATR1
	SaveEMEM(03066, 031102);
	//DAPDATR2
	SaveEMEM(03067, 01111);
	//LEMMASS
	iTemp = SingleToBuffer(BLOCKII.LMMass*LBS2KG, 16);
	SaveEMEM(03073, iTemp);
	//CSMMASS
	iTemp = SingleToBuffer(BLOCKII.CSMMass*LBS2KG, 16);
	SaveEMEM(03074, iTemp);
	//POLYNUM
	SavePOLYNUM(03261);
	//SATRLRT
	SaveEMEM(03300, 0);
	SaveEMEM(03301, 016441);
	//RPSTART
	dTemp = BLOCKII.RPSTART*100.0;
	iTemp = SingleToBuffer(dTemp, 14);
	SaveEMEM(03302, iTemp);
	//POLYSTOP
	dTemp = -BLOCKII.POLYSTOP*100.0;
	iTemp = SingleToBuffer(dTemp, 14);
	SaveEMEM(03303, iTemp);
	//SATRATE
	SaveEMEM(03321, 0);
	SaveEMEM(03322, 0344);
	SaveEMEM(03323, 077433);
	SaveEMEM(03324, 0);
	//SATSCALE
	SaveEMEM(03331, 010000); //0.3 V/DEG
	//P37RANGE
	dTemp = BLOCKII.P37RANGE / 21600.0;
	iTemp = SingleToBuffer(dTemp, 0);
	SaveEMEM(03376, iTemp);
	//ECSTEER
	SaveEMEM(03424, 010000);

	//TB6JOB
	SaveEMEM(03513, 04);
	SaveEMEM(03514, 06);
	SaveEMEM(03515, 031413);
	SaveEMEM(03516, 053140);
	SaveEMEM(03517, 06);
	SaveEMEM(03520, 040025);
	SaveEMEM(03521, 021140);
	SaveEMEM(03522, 05357);
	SaveEMEM(03523, 01550);
	SaveEMEM(03524, 010067);
	SaveEMEM(03525, 011573);
	SaveEMEM(03526, 04106);

	SaveEMEM(03527, 035055);
	SaveEMEM(03530, 05251);
	SaveEMEM(03531, 01543);
	SaveEMEM(03532, 06);
	SaveEMEM(03533, 041413);
	SaveEMEM(03534, 053046);
	SaveEMEM(03535, 01574);
	SaveEMEM(03537, 021046);
	SaveEMEM(03540, 033300);
	SaveEMEM(03541, 04676);
	SaveEMEM(03542, 020707);
	SaveEMEM(03543, 037663);
	SaveEMEM(03544, 05150);

	SaveEMEM(03545, 01525);
	SaveEMEM(03546, 010067);
	SaveEMEM(03547, 05340);
	SaveEMEM(03550, 035017);
	SaveEMEM(03551, 06);
	SaveEMEM(03552, 05012);
	SaveEMEM(03553, 035031);
	SaveEMEM(03554, 06);
	SaveEMEM(03555, 05011);
	SaveEMEM(03556, 06);
	SaveEMEM(03557, 030025);
	SaveEMEM(03560, 053337);
	SaveEMEM(03561, 05303);

	SaveEMEM(03562, 01750);
	SaveEMEM(03563, 045017);
	SaveEMEM(03564, 06);
	SaveEMEM(03565, 03012);
	SaveEMEM(03566, 045031);
	SaveEMEM(03567, 06);
	SaveEMEM(03570, 03011);
	SaveEMEM(03571, 025573);
	SaveEMEM(03572, 05340);
	SaveEMEM(03573, 0);
	SaveEMEM(03574, 06);
	SaveEMEM(03575, 030025);
	SaveEMEM(03576, 01537);
}

void AGCPadloadGenerator::Comanche108Padload()
{
	CMCDefaults(false, true);

	//FLAGWRD1
	SaveEMEM(075, 0);
	//FLAGWRD3
	SaveEMEM(077, 0);
	//FLAGWRD8
	//Bit 8 set: SURFFLAG. LM on lunar surface
	SaveEMEM(0104, 0200);
	//FLAGWRD9
	//Bit 3 set: MID1FLAG. Integrate to TDEC
	SaveEMEM(0105, 04);
	//NO.PASS. Number of passes of P24 before landmark coordinate update.
	SaveEMEM(0736, 037777);
	//PIPTIME
	double PIPTIME = (LaunchMJD - PrelaunchMJD)*8.64e6;
	DoubleToBuffer(PIPTIME, 28, iTemp, iTemp2);
	SaveEMEM(01041, iTemp);
	SaveEMEM(01042, iTemp2);
	//PGNCSALT
	dTemp = PadAlt;
	DoubleToBuffer(dTemp, 29, iTemp, iTemp2);
	SaveEMEM(01133, iTemp);
	SaveEMEM(01134, iTemp2);
	//PADLONG
	dTemp = PadLong / 360.0;
	DoubleToBuffer(dTemp, 0, iTemp, iTemp2);
	SaveEMEM(01135, iTemp);
	SaveEMEM(01136, iTemp2);
	//RTED1. First coefficient defining high speed V-gamma target line polynomial
	DoubleToBuffer(RTED1, 3, iTemp, iTemp2);
	SaveEMEM(01347, iTemp);
	SaveEMEM(01350, iTemp2);
	//DVTHRESH
	iTemp = SingleToBuffer(DVTHRESH / 100.0, -2);
	SaveEMEM(01351, iTemp);
	//Horizalt
	DoubleToBuffer(HORIZALT, 29, iTemp, iTemp2);
	SaveEMEM(01352, iTemp);
	SaveEMEM(01353, iTemp2);
	//ALTVAR
	iTemp = SingleToBuffer(ALTVAR, -16);
	SaveEMEM(01354, iTemp);
	//EMDOT
	dTemp = CMCDATA.EMDOT*LBS2KG;
	iTemp = SingleToBuffer(dTemp / 100.0, 3);
	SaveEMEM(01355, iTemp);
	//EK1VAL
	DoubleToBuffer(CMCDATA.EK1VAL*LBF2N / 100.0, 23, iTemp, iTemp2);
	SaveEMEM(01765, iTemp);
	SaveEMEM(01766, iTemp2);
	//EK2VAL
	DoubleToBuffer(CMCDATA.EK2VAL*LBF2N / 100.0, 23, iTemp, iTemp2);
	SaveEMEM(01767, iTemp);
	SaveEMEM(01770, iTemp2);
	//EK3VAL
	iTemp = SingleToBuffer(CMCDATA.EK3VAL*LBF2N / pow(100.0, 2), 9);
	SaveEMEM(01771, iTemp);
	//FANG
	DoubleToBuffer(CMCDATA.FANG*LBF2N / pow(100.0, 2), 7, iTemp, iTemp2);
	SaveEMEM(01772, iTemp);
	SaveEMEM(01773, iTemp2);

	R2Model(01774);

	//TRUNSF
	iTemp = SingleToBuffer(CMCDATA.TRUNSF, 27);
	SaveEMEM(01776, iTemp);
	//SHAFTSF
	iTemp = SingleToBuffer(CMCDATA.SHAFTSF, 25);
	SaveEMEM(01777, iTemp);

	//LAUNCHAZ
	dTemp = LaunchAzimuth / 360.0;
	DoubleToBuffer(dTemp, 0, iTemp, iTemp2);
	SaveEMEM(02633, iTemp);
	SaveEMEM(02634, iTemp2);

	//WMIDPOS
	dTemp = CMCDATA.WMIDPOS*FT2M;
	iTemp = SingleToBuffer(dTemp, 19);
	SaveEMEM(03000, iTemp);
	//WMIDVEL
	dTemp = CMCDATA.WMIDVEL*FT2M / 100.0;
	iTemp = SingleToBuffer(dTemp, 0);
	SaveEMEM(03001, iTemp);

	//LADPAD
	dTemp = 0.27;
	iTemp = SingleToBuffer(dTemp, 0);
	SaveEMEM(03007, iTemp);
	//LODPAD
	dTemp = 0.207;
	iTemp = SingleToBuffer(dTemp, 0);
	SaveEMEM(03010, iTemp);
	//ALFAPAD
	dTemp = -19.06 / 360.0;
	iTemp = SingleToBuffer(dTemp, -1);
	SaveEMEM(03011, iTemp);
	//P37RANGE
	dTemp = BLOCKII.P37RANGE / 21600.0;
	iTemp = SingleToBuffer(dTemp, 0);
	SaveEMEM(03012, iTemp);
	//ETDECAY
	dTemp = 0.6*100.0;
	iTemp = SingleToBuffer(dTemp, 14);
	SaveEMEM(03013, iTemp);
	//EKPRIME
	SaveEMEM(03014, 0123);
	SaveEMEM(03015, 0175);
	//EKTLX
	SaveEMEM(03016, 017433);
	SaveEMEM(03017, 04500);
	SaveEMEM(03020, 0334);
	//EREPFRAC
	SaveEMEM(03021, 01000);
	SaveEMEM(03022, 0232);
	//PACTOFF
	iTemp = SingleToBuffer(BLOCKII.PACTOFF / 0.023725, 14);
	SaveEMEM(03023, iTemp);
	//YACTOFF
	iTemp = SingleToBuffer(BLOCKII.YACTOFF / 0.023725, 14);
	SaveEMEM(03024, iTemp);
	//HBN10
	SaveEMEM(03025, 037777);
	//HBN11/2
	SaveEMEM(03026, 0);
	//HBN12
	SaveEMEM(03027, 0);
	//HBD11/2
	SaveEMEM(03030, 054360);
	//HBD12
	SaveEMEM(03031, 021075);
	//HBN20
	SaveEMEM(03032, 037777);
	//HBN21/2
	SaveEMEM(03033, 060465);
	//HBN22
	SaveEMEM(03034, 0);
	//HBD21/2
	SaveEMEM(03035, 054360);
	//HBD22
	SaveEMEM(03036, 021075);
	//HBN30
	SaveEMEM(03037, 037777);
	//HBN31/2
	SaveEMEM(03040, 057142);
	//HBN32
	SaveEMEM(03041, 033106);
	//HBD31/2
	SaveEMEM(03042, 050741);
	//HBD32
	SaveEMEM(03043, 031162);
	//DAPDATR1
	SaveEMEM(03065, 031102);
	//DAPDATR2
	SaveEMEM(03066, 01111);
	//LEMMASS
	iTemp = SingleToBuffer(BLOCKII.LMMass*LBS2KG, 16);
	SaveEMEM(03073, iTemp);
	//CSMMASS
	iTemp = SingleToBuffer(BLOCKII.CSMMass*LBS2KG, 16);
	SaveEMEM(03074, iTemp);
	//POLYNUM
	SavePOLYNUM(03261);
	//SATRLRT
	SaveEMEM(03300, 0);
	SaveEMEM(03301, 016441);
	//RPSTART
	dTemp = BLOCKII.RPSTART*100.0;
	iTemp = SingleToBuffer(dTemp, 14);
	SaveEMEM(03302, iTemp);
	//POLYSTOP
	dTemp = -BLOCKII.POLYSTOP*100.0;
	iTemp = SingleToBuffer(dTemp, 14);
	SaveEMEM(03303, iTemp);
	//SATRATE
	SaveEMEM(03321, 0);
	SaveEMEM(03322, 0344);
	SaveEMEM(03323, 077433);
	SaveEMEM(03324, 0);
	//SATSCALE
	SaveEMEM(03331, 010000); //0.3 V/DEG
	//HORISLP
	SaveEMEM(03376, 0);
	SaveEMEM(03377, 0);
}

void AGCPadloadGenerator::Artemis72Padload()
{
	CMCDefaults(false);

	//FLAGWRD1
	SaveEMEM(075, 0);
	//FLAGWRD3
	//Bit 15 set: V50N18FL. Enable R60 att maneuver
	SaveEMEM(077, 040000);
	//FLAGWRD8
	//Bit 8 set: SURFFLAG. LM on lunar surface
	SaveEMEM(0104, 0200);
	//FLAGWRD9
	//Bit 3 set: MID1FLAG. Integrate to TDEC
	SaveEMEM(0105, 04);
	//FLAGWRD10
	SaveEMEM(0106, 0);
	//C31FLWRD
	SaveEMEM(0374, 0);
	//NO.PASS. Number of passes of P24 before landmark coordinate update.
	SaveEMEM(0737, 037777);
	//N26/PRI
	SaveEMEM(01016, 0);
	//N26/2CAD
	SaveEMEM(01017, 0);
	SaveEMEM(01020, 0);
	//PIPTIME
	double PIPTIME = (LaunchMJD - PrelaunchMJD)*8.64e6;
	DoubleToBuffer(PIPTIME, 28, iTemp, iTemp2);
	SaveEMEM(01043, iTemp);
	SaveEMEM(01044, iTemp2);
	//PGNCSALT
	dTemp = PadAlt;
	DoubleToBuffer(dTemp, 29, iTemp, iTemp2);
	SaveEMEM(01133, iTemp);
	SaveEMEM(01134, iTemp2);
	//PADLONG
	dTemp = PadLong / 360.0;
	DoubleToBuffer(dTemp, 0, iTemp, iTemp2);
	SaveEMEM(01135, iTemp);
	SaveEMEM(01136, iTemp2);
	//RTED1. First coefficient defining high speed V-gamma target line polynomial
	DoubleToBuffer(RTED1, 3, iTemp, iTemp2);
	SaveEMEM(01342, iTemp);
	SaveEMEM(01343, iTemp2);
	//DVTHRESH
	iTemp = SingleToBuffer(DVTHRESH / 100.0, -2);
	SaveEMEM(01344, iTemp);
	//Horizalt
	DoubleToBuffer(HORIZALT, 29, iTemp, iTemp2);
	SaveEMEM(01345, iTemp);
	SaveEMEM(01346, iTemp2);
	//ALTVAR
	iTemp = SingleToBuffer(ALTVAR, -16);
	SaveEMEM(01347, iTemp);
	//EMDOT
	dTemp = CMCDATA.EMDOT*LBS2KG;
	iTemp = SingleToBuffer(dTemp / 100.0, 3);
	SaveEMEM(01350, iTemp);
	//EIMP1SEC
	dTemp = CMCDATA.EIMP1SEC*LBF2N / 100.0;
	iTemp = SingleToBuffer(dTemp, 14);
	SaveEMEM(01763, iTemp);
	//EFIMP01
	dTemp = CMCDATA.EFIMP01*LBF2N / 100.0;
	iTemp = SingleToBuffer(dTemp, 14);
	SaveEMEM(01764, iTemp);
	//EFIMP16
	dTemp = CMCDATA.EFIMP16*LBF2N / 100.0;
	iTemp = SingleToBuffer(dTemp, 14);
	SaveEMEM(01765, iTemp);

	R2Model(01766);

	//TRUNSF
	iTemp = SingleToBuffer(CMCDATA.TRUNSF, 27);
	SaveEMEM(01770, iTemp);
	//SHAFTSF
	iTemp = SingleToBuffer(CMCDATA.SHAFTSF, 25);
	SaveEMEM(01771, iTemp);

	//DTF
	dTemp = 0.3*100.0; //0.3 seconds
	iTemp = SingleToBuffer(dTemp, 14);
	SaveEMEM(02354, iTemp);
	//HAMDELH
	dTemp = 15.0*1852.0; //15NM
	DoubleToBuffer(dTemp, 29, iTemp, iTemp2);
	SaveEMEM(02367, iTemp);
	SaveEMEM(02370, iTemp2);
	//WRDTIME
	dTemp = 2400.0*100.0; //2400 seconds
	iTemp = SingleToBuffer(dTemp, 28);
	SaveEMEM(02371, iTemp);
	//MINBLKTM
	dTemp = 328.8*100.0; //328.8 seconds
	iTemp = SingleToBuffer(dTemp, 28);
	SaveEMEM(02372, iTemp);
	//TBEFCCMP
	dTemp = 822.0*100.0; //822 seconds
	iTemp = SingleToBuffer(dTemp, 28);
	SaveEMEM(02373, iTemp);
	//BRNBLKTM
	dTemp = 822.0*100.0; //822 seconds
	iTemp = SingleToBuffer(dTemp, 28);
	SaveEMEM(02374, iTemp);
	//MAXWTIME
	dTemp = 3616.8*100.0; //3616.8 seconds
	iTemp = SingleToBuffer(dTemp, 28);
	SaveEMEM(02375, iTemp);
	//FINCMPTM
	dTemp = 493.2*100.0; //493.2 seconds
	iTemp = SingleToBuffer(dTemp, 28);
	SaveEMEM(02376, iTemp);

	//LAUNCHAZ
	dTemp = LaunchAzimuth / 360.0;
	DoubleToBuffer(dTemp, 0, iTemp, iTemp2);
	SaveEMEM(02633, iTemp);
	SaveEMEM(02634, iTemp2);

	//WMIDPOS
	dTemp = CMCDATA.WMIDPOS*FT2M;
	iTemp = SingleToBuffer(dTemp, 19);
	SaveEMEM(03000, iTemp);
	//WMIDVEL
	dTemp = CMCDATA.WMIDVEL*FT2M / 100.0;
	iTemp = SingleToBuffer(dTemp, 0);
	SaveEMEM(03001, iTemp);

	//LADPAD
	dTemp = BLOCKII.LADPAD;
	iTemp = SingleToBuffer(dTemp, 0);
	SaveEMEM(03007, iTemp);
	//LODPAD
	dTemp = BLOCKII.LODPAD;
	iTemp = SingleToBuffer(dTemp, 0);
	SaveEMEM(03010, iTemp);
	//ALFAPAD
	dTemp = BLOCKII.ALFAPAD / 360.0;
	iTemp = SingleToBuffer(dTemp, -1);
	SaveEMEM(03011, iTemp);
	//P37RANGE
	dTemp = BLOCKII.P37RANGE / 21600.0;
	iTemp = SingleToBuffer(dTemp, 0);
	SaveEMEM(03012, iTemp);
	//ETDECAY
	dTemp = 0.6*100.0;
	iTemp = SingleToBuffer(dTemp, 14);
	SaveEMEM(03013, iTemp);
	//EKPRIME
	SaveEMEM(03014, 0123);
	SaveEMEM(03015, 0175);
	//EKTLX
	SaveEMEM(03016, 017433);
	SaveEMEM(03017, 04500);
	SaveEMEM(03020, 0334);
	//EREPFRAC
	SaveEMEM(03021, 01000);
	SaveEMEM(03022, 0232);
	//PACTOFF
	iTemp = SingleToBuffer(BLOCKII.PACTOFF / 0.023725, 14);
	SaveEMEM(03023, iTemp);
	//YACTOFF
	iTemp = SingleToBuffer(BLOCKII.YACTOFF / 0.023725, 14);
	SaveEMEM(03024, iTemp);
	//HBN10
	SaveEMEM(03025, 037777);
	//HBN11/2
	SaveEMEM(03026, 0);
	//HBN12
	SaveEMEM(03027, 0);
	//HBD11/2
	SaveEMEM(03030, 054360);
	//HBD12
	SaveEMEM(03031, 021075);
	//HBN20
	SaveEMEM(03032, 037777);
	//HBN21/2
	SaveEMEM(03033, 060465);
	//HBN22
	SaveEMEM(03034, 0);
	//HBD21/2
	SaveEMEM(03035, 054360);
	//HBD22
	SaveEMEM(03036, 021075);
	//HBN30
	SaveEMEM(03037, 037777);
	//HBN31/2
	SaveEMEM(03040, 057142);
	//HBN32
	SaveEMEM(03041, 033106);
	//HBD31/2
	SaveEMEM(03042, 050741);
	//HBD32
	SaveEMEM(03043, 031162);
	//DAPDATR1
	SaveEMEM(03065, 031102);
	//DAPDATR2
	SaveEMEM(03066, 01111);
	//LEMMASS
	iTemp = SingleToBuffer(BLOCKII.LMMass*LBS2KG, 16);
	SaveEMEM(03072, iTemp);
	//CSMMASS
	iTemp = SingleToBuffer(BLOCKII.CSMMass*LBS2KG, 16);
	SaveEMEM(03073, iTemp);
	//POLYNUM
	SavePOLYNUM(03261);
	//SATRLRT
	SaveEMEM(03300, 0);
	SaveEMEM(03301, 016441);
	//RPSTART
	dTemp = BLOCKII.RPSTART*100.0;
	iTemp = SingleToBuffer(dTemp, 14);
	SaveEMEM(03302, iTemp);
	//POLYSTOP
	dTemp = -BLOCKII.POLYSTOP*100.0;
	iTemp = SingleToBuffer(dTemp, 14);
	SaveEMEM(03303, iTemp);
	//SATRATE
	SaveEMEM(03310, 0);
	SaveEMEM(03311, 0344);
	SaveEMEM(03312, 077433);
	SaveEMEM(03313, 0);
	//SATSCALE
	SaveEMEM(03320, 010000); //0.3 V/DEG
	//HORISLP
	SaveEMEM(03375, 0);
	SaveEMEM(03376, 0);
}

void AGCPadloadGenerator::SavePOLYNUM(int address)
{
	//POLYNUM
	SaveEMEM(address++, 05);

	DoubleToBuffer(BLOCKII.POLYNUM[0] / 360.0 / pow(100.0, 0), 5, iTemp, iTemp2);
	SaveEMEM(address++, iTemp);
	SaveEMEM(address++, iTemp2);

	DoubleToBuffer(BLOCKII.POLYNUM[1] / 360.0 / pow(100.0, 1), -9, iTemp, iTemp2);
	SaveEMEM(address++, iTemp);
	SaveEMEM(address++, iTemp2);

	DoubleToBuffer(BLOCKII.POLYNUM[2] / 360.0 / pow(100.0, 2), -23, iTemp, iTemp2);
	SaveEMEM(address++, iTemp);
	SaveEMEM(address++, iTemp2);

	DoubleToBuffer(BLOCKII.POLYNUM[3] / 360.0 / pow(100.0, 3), -37, iTemp, iTemp2);
	SaveEMEM(address++, iTemp);
	SaveEMEM(address++, iTemp2);

	DoubleToBuffer(BLOCKII.POLYNUM[4] / 360.0 / pow(100.0, 4), -51, iTemp, iTemp2);
	SaveEMEM(address++, iTemp);
	SaveEMEM(address++, iTemp2);

	DoubleToBuffer(BLOCKII.POLYNUM[5] / 360.0 / pow(100.0, 5), -65, iTemp, iTemp2);
	SaveEMEM(address++, iTemp);
	SaveEMEM(address++, iTemp2);

	DoubleToBuffer(BLOCKII.POLYNUM[6] / 360.0 / pow(100.0, 6), -79, iTemp, iTemp2);
	SaveEMEM(address++, iTemp);
	SaveEMEM(address++, iTemp2);
}

MATRIX3 AGCPadloadGenerator::SolariumEarthFixedToSM(double lat, double lng, double azi)
{
	VECTOR3 R_P = r_from_latlong(lat, lng);
	VECTOR3 REFS6 = unit(-R_P);
	VECTOR3 E = unit(crossp(REFS6, _V(0, 0, 1)));
	VECTOR3 S = unit(crossp(E, REFS6));
	VECTOR3 REFS0 = E * sin(azi) + S * cos(azi);
	VECTOR3 REFS3 = unit(crossp(REFS6, REFS0));
	MATRIX3 REFS = _M(-REFS6.x, -REFS6.y, -REFS6.z, REFS3.x, REFS3.y, REFS3.z, REFS0.x, REFS0.y, REFS0.z);
	return REFS;
}

double AGCPadloadGenerator::Solarium055DTEPOCHCalculation(double A_Z0, double MJD_0, double MJD_L, double lng)
{
	double t0 = (MJD_L - MJD_0) * 24 * 3600;
	return fmod(A_Z0 + w_Earth * t0 + lng, PI2) / w_Earth;
}

double AGCPadloadGenerator::HANGLE(int E, int Y, int D)
{
	int XN;
	static const double A = 0.0929;
	static const double B = 8640184.542;
	static const double W1 = 1.720217954160054e-2;
	double C, T, DE, BHA, DI, DELTA;

	XN = (E - 1901) / 4;
	C = -86400.0*(double)(E - 1900) - 74.164;
	T = 2 * C / (-B - sqrt(B*B - 4 * A*C));
	DE = 36525.0*T - 365.0*(double)(E - 1900) + 0.5 - (double)XN;
	if (Y == E)
	{
		DI = D;
	}
	else
	{
		int X = Y % 4;
		if (X == 0)
		{
			DI = D - 366.0;
		}
		else
		{
			DI = D - 365.0;
		}
	}
	DELTA = DI - DE;
	BHA = PI2 / 3.6 + W1 * DELTA;
	return BHA;
}

MATRIX3 AGCPadloadGenerator::CalculateEarthTransformationMatrix(double t_M, double A_Z0, double w_E)
{
	double A_Z; 

	A_Z = A_Z0 + w_E * t_M;

	return _M(cos(A_Z), sin(A_Z), 0., -sin(A_Z), cos(A_Z), 0., 0., 0., 1.);
}

MATRIX3 AGCPadloadGenerator::CalculateMoonTransformationMatrix(double t_M, double B_0, double B_dot, double Omega_I0, double Omega_I_dot, double F_0, double F_dot, double cosI, double sinI)
{
	MATRIX3 M1, M2, M3, M4;
	double B, Omega_I, F;

	B = B_0 + B_dot * t_M;
	Omega_I = Omega_I0 + Omega_I_dot * t_M;
	F = F_0 + F_dot * t_M;
	M1 = _M(1., 0., 0., 0., cos(B), sin(B), 0., -sin(B), cos(B));
	M2 = _M(cos(Omega_I), sin(Omega_I), 0., -sin(Omega_I), cos(Omega_I), 0., 0., 0., 1.);
	M3 = _M(1., 0., 0., 0., cosI, -sinI, 0., sinI, cosI);
	M4 = _M(-cos(F), -sin(F), 0., sin(F), -cos(F), 0., 0., 0., 1.);
	return mul(M4, mul(M3, mul(M2, M1)));
}

AGCPadloadGenerator::AGCVersions AGCPadloadGenerator::GetCMCVersion(std::string name)
{
	if (name == "Colossus237")
	{
		return AGCVersions::Colossus237;
	}
	else if (name == "Colossus249")
	{
		return AGCVersions::Colossus249;
	}
	else if (name == "Comanche045")
	{
		return AGCVersions::Comanche045;
	}
	else if (name == "Comanche055")
	{
		return AGCVersions::Comanche055;
	}
	else if (name == "Comanche067")
	{
		return AGCVersions::Comanche067;
	}
	else if (name == "Comanche072")
	{
		return AGCVersions::Comanche072;
	}
	else if (name == "Artemis072")
	{
		return AGCVersions::Artemis072;
	}
	else if (name == "Skylark048")
	{
		return AGCVersions::Skylark048;
	}

	return AGCVersions::AGCVersionError;
}

AGCPadloadGenerator::AGCVersions AGCPadloadGenerator::GetLGCVersion(std::string name)
{
	if (name == "Sundance306")
	{
		return AGCVersions::Sundance306;
	}
	else if (name == "Luminary069")
	{
		return AGCVersions::Luminary069;
	}
	else if (name == "Luminary069R2")
	{
		return AGCVersions::Luminary069R2;
	}
	else if (name == "Luminary099")
	{
		return AGCVersions::Luminary099;
	}
	else if (name == "Luminary116")
	{
		return AGCVersions::Luminary116;
	}
	else if (name == "Luminary131")
	{
		return AGCVersions::Luminary131;
	}
	else if (name == "Luminary131R1")
	{
		return AGCVersions::Luminary131R1;
	}
	else if (name == "Luminary178")
	{
		return AGCVersions::Luminary178;
	}
	else if (name == "Zerlina56")
	{
		return AGCVersions::Zerlina56;
	}
	else if (name == "Luminary210")
	{
		return AGCVersions::Luminary210;
	}

	return AGCVersions::AGCVersionError;
}

int AGCPadloadGenerator::GetPIOSDataSet(std::string name, PIOSDataSet &data)
{
	std::ifstream file;
	std::string line, linename;
	char Buffer[256];
	int inttemp[2];

	file.open("PIOSDataSets.txt");

	if (file.is_open() == false) return 1;

	while (std::getline(file, line))
	{
		if (sscanf_s(line.c_str(), "%s", Buffer, 255) == 1)
		{
			linename.assign(Buffer);

			if (name == linename)
			{
				//Load
				if (sscanf_s(line.c_str(), "%s %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf %d", Buffer, 255, &data.epoch, &data.w_E, &data.B_0, &data.Omega_I0, &data.F_0, &data.B_dot, &data.Omega_I_dot,
					&data.F_dot, &data.cosI, &data.sinI, &data.t0, &inttemp[0], &data.AZ0, &inttemp[1]) != 15)
				{
					return 1;
				}

				data.name.assign(Buffer);
				data.AZ0Hardcoded = (inttemp[0] != 0);
				data.ExtendedLimit = (inttemp[1] != 0);
				break;
			}
		}
	}

	file.close();

	return 0;
}

void AGCPadloadGenerator::AGCCorrectionVectors(PIOSDataSet dataset, double mjd_launchday, double dt_UNITW, double dt_504LM, bool IsCMC, bool IsSundance)
{
	MATRIX3 R, Rot2, J2000, R2, R3, M, M_AGC;
	VECTOR3 UNITW;
	double mjd_UNITW, brcsmjd;
	double A_Z, A_Z0, mjd_504LM;
	int mem;

	mjd_UNITW = mjd_launchday + dt_UNITW;
	mjd_504LM = mjd_launchday + dt_504LM;

	Rot2 = _M(1., 0., 0., 0., 0., 1., 0., 1., 0.);
	R = GetRotationMatrix(BODY_EARTH, mjd_UNITW);

	//EARTH ROTATIONS
	brcsmjd = MJDOfNBYEpoch(dataset.epoch);
	J2000 = J2000EclToBRCSMJD(brcsmjd);
	R2 = mul(tmat(Rot2), mul(R, Rot2));
	R3 = mul(J2000, R2);

	UNITW = mul(R3, _V(0, 0, 1));

	A_Z = atan2(R3.m21, R3.m11);

	//Hardcoded A_Z0?
	if (dataset.AZ0Hardcoded)
	{
		A_Z0 = dataset.AZ0;
	}
	else
	{
		A_Z0 = fmod((A_Z - dataset.w_E * (mjd_UNITW - dataset.t0) * 24.0 * 3600.0), PI2);  //AZ0 for mission
		if (A_Z0 < 0) A_Z0 += PI2;
	}

	M = _M(cos(A_Z), sin(A_Z), 0., -sin(A_Z), cos(A_Z), 0., 0., 0., 1.);
	M_AGC = mul(R3, M);

	mem = 01711;

	if (dataset.AZ0Hardcoded == false)
	{
		DoubleToBuffer(A_Z0 / PI2, 0, iTemp, iTemp2);

		SaveEMEM(mem, iTemp); mem++;
		SaveEMEM(mem, iTemp2); mem++;
	}

	DoubleToBuffer(UNITW.x, 0, iTemp, iTemp2);
	SaveEMEM(mem, iTemp); mem++;
	SaveEMEM(mem, iTemp2); mem++;

	DoubleToBuffer(UNITW.y, 0, iTemp, iTemp2);
	SaveEMEM(mem, iTemp); mem++;
	SaveEMEM(mem, iTemp2); mem++;
	if (IsCMC || IsSundance)
	{
		DoubleToBuffer(UNITW.z, 0, iTemp, iTemp2);
		SaveEMEM(mem, iTemp); mem++;
		SaveEMEM(mem, iTemp2); mem++;
	}

	//UNITW deviation analysis

	debugfile.open("DebugData.txt");
	debugfile << "CORRECTION VECTORS\n\n";
	debugfile << "Earth correction vector deviation analysis\n\n";

	MATRIX3 RM, MM;
	VECTOR3 RLS_Earth,  RLS_ecl, RLS_BRCS, RLS_BRCS2, l, l_P;
	double MJD, err, t_M;

	MJD = mjd_UNITW;
	l = _V(-UNITW.y, UNITW.x, 0.0);
	RLS_Earth = r_from_latlong(28.5*RAD, -80.5*RAD, 0.0, 0, 0);

	while (MJD < mjd_UNITW + 14.5)
	{
		t_M = (MJD - dataset.t0) * 24.0 * 3600.0;
		RM = GetRotationMatrix(BODY_EARTH, MJD);

		//Orbiter
		RLS_ecl = rhmul(RM, RLS_Earth);
		RLS_BRCS = mul(J2000, RLS_ecl);

		//AGC
		MM = CalculateEarthTransformationMatrix(t_M, A_Z0, dataset.w_E);
		l_P = mul(MM, l);
		RLS_BRCS2 = tmul(MM, RLS_Earth + crossp(l_P, RLS_Earth));

		err = acos(dotp(unit(RLS_BRCS), unit(RLS_BRCS2)))*R_Earth;

		sprintf_s(Buffer, 256, "%.3f = %g (meters)\n", MJD, err);
		debugfile << Buffer;

		MJD += 0.5;
	}

	if (IsSundance)
	{
		debugfile << std::endl;
		debugfile.close();
		return;
	}

	//MOON ROTATIONS
	MATRIX3 M_AGC_M, R2M, R3M;
	VECTOR3 lm;

	t_M = (mjd_504LM - dataset.t0) * 24.0 * 3600.0;
	RM = GetRotationMatrix(BODY_MOON, mjd_504LM);
	MM = CalculateMoonTransformationMatrix(t_M, dataset.B_0, dataset.B_dot, dataset.Omega_I0, dataset.Omega_I_dot, dataset.F_0, dataset.F_dot, dataset.cosI, dataset.sinI);

	R2M = mul(tmat(Rot2), mul(RM, Rot2));
	R3M = mul(J2000, R2M);
	M_AGC_M = mul(MM, R3M);

	lm.x = atan2(M_AGC_M.m32, M_AGC_M.m33);
	lm.y = atan2(-M_AGC_M.m31, sqrt(M_AGC_M.m32*M_AGC_M.m32 + M_AGC_M.m33*M_AGC_M.m33));
	lm.z = atan2(M_AGC_M.m21, M_AGC_M.m11);

	if (IsCMC)
	{
		mem = 02011;
	}
	else
	{
		mem = 02012;
	}

	DoubleToBuffer(lm.x, 0, iTemp, iTemp2);
	SaveEMEM(mem, iTemp); mem++;
	SaveEMEM(mem, iTemp2); mem++;
	DoubleToBuffer(lm.y, 0, iTemp, iTemp2);
	SaveEMEM(mem, iTemp); mem++;
	SaveEMEM(mem, iTemp2); mem++;
	DoubleToBuffer(lm.z, 0, iTemp, iTemp2);
	SaveEMEM(mem, iTemp); mem++;
	SaveEMEM(mem, iTemp2); mem++;

	debugfile << "\nLunar correction vector deviation analysis\n\n";

	sprintf_s(Buffer, 256, "504LM = %lf %lf %lf\n\n", lm.x, lm.y, lm.z);
	debugfile << Buffer;

	//504LM deviation analysis
	MJD = mjd_504LM - 3.0;

	while (MJD < mjd_504LM + 3.0)
	{
		t_M = (MJD - dataset.t0) * 24.0 * 3600.0;
		RM = GetRotationMatrix(BODY_MOON, MJD);

		//Orbiter
		RLS_ecl = rhmul(RM, RLS);
		RLS_BRCS = mul(J2000, RLS_ecl);

		//AGC
		MM = CalculateMoonTransformationMatrix(t_M, dataset.B_0, dataset.B_dot, dataset.Omega_I0, dataset.Omega_I_dot, dataset.F_0, dataset.F_dot, dataset.cosI, dataset.sinI);
		RLS_BRCS2 = tmul(MM, RLS + crossp(lm, RLS));

		err = acos2(dotp(unit(RLS_BRCS), unit(RLS_BRCS2)))*R_Moon;

		sprintf_s(Buffer, 256, "%.3f = %g (meters)\n", MJD, err);
		debugfile << Buffer;

		MJD += 0.5;
	}

	debugfile << std::endl;
	debugfile.close();
}

void AGCPadloadGenerator::AGCEphemeris(double T0, int Epoch, double TEphem0, double Span)
{
	//T0: Start time of ephemeris (MJD)
	//Epoch: Year of BRCS
	//TEphem0: Time reference of yearly coordinate system (MJD)
	//Span: Ephemeris span in days

	const int datapoints = 57;
	const int poly = 10;
	const double JulianDateDT = 2400000.5;

	//Calculate mid point
	double Tm0 = T0 + Span / 2.0;
	//Calculate step length
	double dt = Span / ((double)(datapoints - 1));
	double T_cur = T0;

	double x[datapoints];
	double y[datapoints];
	double xx[poly];
	double yy[poly];
	double zz[poly];

	VECTOR3 pos[datapoints];

	double dscale = 1.0 / pow(2.0, 31.0);
	double tscale = 8.64e6 / pow(2.0, 26.0);

	double brcsmjd = MJDOfNBYEpoch(Epoch);
	MATRIX3 Mat = J2000EclToBRCSMJD(brcsmjd);//J2000EclToEqu(Epoch);

	// ----- Solar Ephemeris -----

	VECTOR3 _SPos1, _SVel;

	agcCelBody_RH(BODY_EARTH, Tm0, EPHEM_TRUEPOS | EPHEM_TRUEVEL, &_SPos1, &_SVel);

	_SPos1 = mul(Mat, -_SPos1);
	_SVel = mul(Mat, -_SVel);

	double avel1 = length(crossp(unit(_SPos1), _SVel / length(_SPos1))) / 100.0;	// Angular velocity in rad/cs
	double avel2 = pow(2.0, 26.0)*avel1 / PI2;										// Convert to AGC scaling

	// ----- Lunar Ephemeris -----

	VECTOR3 _LPos;

	agcCelBody_RH(BODY_MOON, Tm0, EPHEM_TRUEPOS, &_LPos);

	for (int i = 0;i < datapoints;i++) {
		x[i] = (T_cur - Tm0) * tscale;
		VECTOR3 _x;
		agcCelBody_RH(BODY_MOON, T_cur, EPHEM_TRUEPOS, &_x);
		pos[i] = mul(Mat, _x) * dscale;
		T_cur += dt;
	}

	for (int i = 0;i < datapoints;i++) y[i] = pos[i].x;
	SolveSeries(x, y, datapoints, xx, poly);
	for (int i = 0;i < datapoints;i++) y[i] = pos[i].y;
	SolveSeries(x, y, datapoints, yy, poly);
	for (int i = 0;i < datapoints;i++) y[i] = pos[i].z;
	SolveSeries(x, y, datapoints, zz, poly);

	int mem = 02033;

	debugfile.open("DebugData.txt", std::ofstream::out | std::ofstream::app);

	debugfile << "EPHEMERIDES\n\n";

	sprintf_s(Buffer, 256, "Epoch   = %d (Year) Epoch of Basic Reference Coordinate System\n", Epoch);
	debugfile << Buffer;
	sprintf_s(Buffer, 256, "Epoch   = %6.6f (MJD) %6.6f (UJD) Epoch of Basic Reference Coordinate System\n", brcsmjd, brcsmjd + JulianDateDT);
	debugfile << Buffer;
	sprintf_s(Buffer, 256, "TEphem0 = %6.6f (MJD) %6.6f (UJD) Ephemeris Time Zero\n", TEphem0, TEphem0 + JulianDateDT);
	debugfile << Buffer;
	sprintf_s(Buffer, 256, "T0      = %6.6f (MJD) %6.6f (EJD) Ephemeris Start Time\n", T0, T0 + JulianDateDT);
	debugfile << Buffer;
	sprintf_s(Buffer, 256, "TM0     = %6.6f (MJD) %6.6f (EJD) Ephemeris Center Time\n", Tm0, Tm0 + JulianDateDT);
	debugfile << Buffer;


	//fprintf(file, "------- TIMEM0 -------\n");

	double t = (Tm0 - TEphem0) * 8640000.0;

	TripleToBuffer(t, 42, iTemp, iTemp2, iTemp3);
	SaveEMEM(mem, iTemp); mem++;
	SaveEMEM(mem, iTemp2); mem++;
	SaveEMEM(mem, iTemp3); mem++;


	//fprintf(file, "------- Lunar Ephemeris -------\n");

	for (int i = 9;i >= 0;i--) {
		DoubleToBuffer(xx[i], 0, iTemp, iTemp2);
		SaveEMEM(mem, iTemp); mem++;
		SaveEMEM(mem, iTemp2); mem++;
		DoubleToBuffer(yy[i], 0, iTemp, iTemp2);
		SaveEMEM(mem, iTemp); mem++;
		SaveEMEM(mem, iTemp2); mem++;
		DoubleToBuffer(zz[i], 0, iTemp, iTemp2);
		SaveEMEM(mem, iTemp); mem++;
		SaveEMEM(mem, iTemp2); mem++;
	}
	//fprintf(file, "------- Solar Ephemeris -------\n");

	mem = 02132;

	DoubleToBuffer(_SPos1.x, 38, iTemp, iTemp2);
	SaveEMEM(mem, iTemp); mem++;
	SaveEMEM(mem, iTemp2); mem++;
	DoubleToBuffer(_SPos1.y, 38, iTemp, iTemp2);
	SaveEMEM(mem, iTemp); mem++;
	SaveEMEM(mem, iTemp2); mem++;
	DoubleToBuffer(_SPos1.z, 38, iTemp, iTemp2);
	SaveEMEM(mem, iTemp); mem++;
	SaveEMEM(mem, iTemp2); mem++;

	DoubleToBuffer(_SVel.x / 100.0, 9, iTemp, iTemp2);
	SaveEMEM(mem, iTemp); mem++;
	SaveEMEM(mem, iTemp2); mem++;
	DoubleToBuffer(_SVel.y / 100.0, 9, iTemp, iTemp2);
	SaveEMEM(mem, iTemp); mem++;
	SaveEMEM(mem, iTemp2); mem++;
	DoubleToBuffer(_SVel.z / 100.0, 9, iTemp, iTemp2);
	SaveEMEM(mem, iTemp); mem++;
	SaveEMEM(mem, iTemp2); mem++;

	DoubleToBuffer(avel2, 0, iTemp, iTemp2);
	SaveEMEM(mem, iTemp); mem++;
	SaveEMEM(mem, iTemp2); mem++;

	//
	// Lunar Ephem Analysis ----------------------------------------------------
	//

	sprintf_s(Buffer, 256, "\nLunar ephemeris deviation analysis\n\n");
	debugfile << Buffer;

	for (int i = 0;i < datapoints;i++) {
		double a = xx[0]; double b = yy[0];
		double c = zz[0]; double h = x[i];

		for (int k = 1;k < 10;k++) { a += xx[k] * h; b += yy[k] * h; c += zz[k] * h; h *= x[i]; }

		VECTOR3 _g = _V(a, b, c);
		double m = Tm0 + x[i] / tscale;

		sprintf_s(Buffer, 256, "%.3f = %g (arcsec)\n", m, DEG*angle(pos[i], _g)*3600.0);
		debugfile << Buffer;
	}

	debugfile << std::endl << std::endl;

	//
	// Solar Ephem Analysis ----------------------------------------------------
	//

	sprintf_s(Buffer, 256, "\nSolar ephemeris deviation analysis\n\n");
	debugfile << Buffer;

	double poserr;

	for (int i = 0;i < datapoints;i++) {
		VECTOR3 _g, _g2;
		double t_M;
		double m = Tm0 + x[i] / tscale;
		agcCelBody_RH(BODY_EARTH, m, EPHEM_TRUEPOS, &_g);
		_g = mul(Mat, -_g);

		t_M = (m - TEphem0) * 8640000.0 - t;
		_g2 = _SPos1 * cos(avel1*t_M) + crossp(_SPos1, unit(crossp(_SVel, _SPos1)))*sin(avel1*t_M);
		poserr = length(_g - _g2);

		sprintf_s(Buffer, 256, "%.3f = %g (arcsec), %g (kilometers)\n", m, DEG*angle(_g, _g2)*3600.0, poserr / 1000.0);
		debugfile << Buffer;
	}

	sprintf_s(Buffer, 256, "\n------- State Vectors at TIMEM0 -------\n");
	debugfile << Buffer;
	sprintf_s(Buffer, 256, "Moon_x = %6.6f\n", _LPos.x);
	debugfile << Buffer;
	sprintf_s(Buffer, 256, "Moon_y = %6.6f\n", _LPos.y);

	sprintf_s(Buffer, 256, "Moon_z = %6.6f\n", _LPos.z);
	debugfile << Buffer;

	sprintf_s(Buffer, 256, "EMBC_x = %6.6f\n", _SPos1.x);
	debugfile << Buffer;
	sprintf_s(Buffer, 256, "EMBC_y = %6.6f\n", _SPos1.y);
	debugfile << Buffer;
	sprintf_s(Buffer, 256, "EMBC_z = %6.6f\n\n", _SPos1.z);
	debugfile << Buffer;

	debugfile.close();
}

int AGCPadloadGenerator::clbkEphemeris(int body, double mjd, int req, double *ret)
{
	if (body == BODY_MOON)
	{
		ELP82(mjd, ret);
		if (req & (EPHEM_BARYPOS | EPHEM_BARYVEL))
			for (int i = 6; i < 12; i++) ret[i] = ret[i - 6];
		return req | (EPHEM_TRUEPOS | EPHEM_TRUEVEL | EPHEM_BARYISTRUE);
	}
	else
	{
		return earth.clbkEphemeris(mjd, req, ret);
	}
}

int AGCPadloadGenerator::agcCelBody_LH(int Cel, double mjd, int Flags, VECTOR3 *Pos, VECTOR3 *Vel)
{
	//
	// Local Variables and structures
	//
	struct s_ret {
		VECTOR3 TrPos;
		VECTOR3 TrVel;
		VECTOR3 BCPos;
		VECTOR3 BCVel;
	} ret;


	int info = clbkEphemeris(Cel, mjd, Flags, (double *)&ret);
	//
	// Data Conversions
	//
	if (info&EPHEM_POLAR) {
		if (info&EPHEM_TRUEVEL && info&EPHEM_TRUEPOS) agcPolar2Cartesian_LH(&ret.TrPos, &ret.TrVel);
		else if (info&EPHEM_TRUEPOS)				  agcPolar2Cartesian_LH(&ret.TrPos, NULL);

		if (info&EPHEM_BARYVEL && info&EPHEM_BARYPOS) agcPolar2Cartesian_LH(&ret.BCPos, &ret.BCVel);
		else if (info&EPHEM_TRUEPOS)                  agcPolar2Cartesian_LH(&ret.BCPos, NULL);
	}


	if (Flags&EPHEM_TRUEPOS) {

		if (info&EPHEM_TRUEPOS) {
			if (Pos) *Pos = ret.TrPos;
			if (Vel) *Vel = ret.TrVel;
			return info & (EPHEM_TRUEPOS | EPHEM_TRUEVEL);
		}
		else {
			if (Pos) *Pos = ret.BCPos;
			if (Vel) *Vel = ret.BCVel;
			return info & (EPHEM_BARYPOS | EPHEM_BARYVEL);
		}
	}

	if (Flags&EPHEM_BARYPOS) {

		if (info&EPHEM_BARYPOS) {
			if (Pos) *Pos = ret.BCPos;
			if (Vel) *Vel = ret.BCVel;
			return info & (EPHEM_BARYPOS | EPHEM_BARYVEL);
		}
		else {
			if (Pos) *Pos = ret.TrPos;
			if (Vel) *Vel = ret.TrVel;
			return info & (EPHEM_TRUEPOS | EPHEM_TRUEVEL);
		}
	}

	if (info&EPHEM_BARYPOS) {
		if (Pos) *Pos = ret.BCPos;
		if (Vel) *Vel = ret.BCVel;
		return info & (EPHEM_BARYPOS | EPHEM_BARYVEL);
	}

	if (Pos) *Pos = ret.TrPos;
	if (Vel) *Vel = ret.TrVel;
	return info & (EPHEM_TRUEPOS | EPHEM_TRUEVEL);
}

int AGCPadloadGenerator::agcCelBody_RH(int Cel, double mjd, int Flags, VECTOR3 *Pos, VECTOR3 *Vel)
{
	int ret = agcCelBody_LH(Cel, mjd, Flags, Pos, Vel);
	if (Pos) agcSwap(Pos);
	if (Vel) agcSwap(Vel);
	return ret;
}

VECTOR3 AGCPadloadGenerator::r_from_latlong(double lat, double lng)
{
	return _V(cos(lng)*cos(lat), sin(lng)*cos(lat), sin(lat));
}

VECTOR3 AGCPadloadGenerator::r_from_latlong(double lat, double lng, double alt, int P, int F)
{
	//P = 0 for earth, 1 = for moon
	//F = 0 for launch pad or landing site radius, 1 = for Fisher ellipsoid or mean lunar radius

	VECTOR3 R_P;
	double gamma, SINL, r0;

	if (P == BODY_EARTH)
	{
		gamma = pow(b_Fisher / a_Fisher, 2);
	}
	else
	{
		gamma = 1.0;
	}

	R_P = unit(_V(cos(lng)*cos(lat), sin(lng)*cos(lat), gamma*sin(lat)));
	SINL = R_P.z;

	if (P == BODY_EARTH)
	{
		if (F == 1)
		{
			r0 = sqrt(b_Fisher*b_Fisher / (1.0 - (1.0 - pow(b_Fisher / a_Fisher, 2))*(1.0 - pow(SINL, 2))));
		}
		else
		{
			r0 = R_Earth;
		}
	}
	else
	{
		r0 = R_Moon; //TBD: landing site radius
	}
	return R_P * (r0 + alt);
}

void AGCPadloadGenerator::BlockIDefaults()
{
	//IMU COMPENSATION
	//Gyro bias drift
	//GBIASX
	SaveEMEM(0744, 0);
	//GBIASY
	SaveEMEM(0745, 0);
	//GBIASZ
	SaveEMEM(0746, 0);
	//Acceleration-sensitive gyro drift along the input axis
	//ADIAX
	SaveEMEM(0747, 0);
	//ADIAY
	SaveEMEM(0750, 0);
	//ADIAZ
	SaveEMEM(0751, 0);
	//Acceleration-sensitive gyro drift along the sin-reference axis
	//ADSRAX
	SaveEMEM(0752, 0);
	//ADSRAY
	SaveEMEM(0753, 0);
	//ADSRAZ
	SaveEMEM(0754, 0);
	//PIPA bias factor
	//PBIASX
	SaveEMEM(0736, 0);
	//PBIASY
	SaveEMEM(0740, 0);
	//PBIASZ
	SaveEMEM(0742, 0);
	//PIPA scale factor
	//PIPASCFX
	SaveEMEM(0737, 0);
	//PIPASCFY
	SaveEMEM(0741, 0);
	//PIPASCFZ
	SaveEMEM(0743, 0);

	//PRELAUNCH ALIGNMENT
	//DTEPOCH
	DoubleToBuffer(BLOCKI.DTEPOCH*100.0, 28, iTemp, iTemp2);
	SaveEMEM(01073, iTemp);
	SaveEMEM(01074, iTemp2);
	//VAZ - Azimuth of vehicle z-axis east of north
	dTemp = PadAzimuth;
	DoubleToBuffer(dTemp / 360.0, 0, iTemp, iTemp2);
	SaveEMEM(01352, iTemp);
	SaveEMEM(01353, iTemp2);
	//LATITUDE - Geodetic latitude of launch pad
	DoubleToBuffer(PadLat / 360.0, 0, iTemp, iTemp2);
	SaveEMEM(01314, iTemp);
	SaveEMEM(01315, iTemp2);
	//AZIMUTH - Azimuth of ??? z-axis east of north
	DoubleToBuffer(LaunchAzimuth / 360.0, 0, iTemp, iTemp2);
	SaveEMEM(01316, iTemp);
	SaveEMEM(01317, iTemp2);
	//POLYENTR - Transfer to ??? routine
	SaveEMEM(01573, 05554);
	//POLYEND - Transfer at end of polynomial
	SaveEMEM(01613, 04024);

	//TROLL - Time from liftoff at which roll monitor begins
	dTemp = BLOCKI.TROLL; //Seconds
	DoubleToBuffer(dTemp*100.0, 28, iTemp, iTemp2);
	SaveEMEM(01562, iTemp);
	SaveEMEM(01563, iTemp2);

	//TPITCH - Time from liftoff at which pitch monitor begins
	dTemp = BLOCKI.TPITCH; //Seconds
	DoubleToBuffer(dTemp*100.0, 28, iTemp, iTemp2);
	SaveEMEM(01564, iTemp);
	SaveEMEM(01565, iTemp2);

	//TENDPITCH - Time pitch monitor is on
	dTemp = BLOCKI.TENDPITCH; //Seconds
	DoubleToBuffer(dTemp*100.0, 28, iTemp, iTemp2);
	SaveEMEM(01566, iTemp);
	SaveEMEM(01567, iTemp2);

	//TTUMON - Nominal time from end of pitch monitor to start of tumble monitor
	dTemp = 43.35; //Seconds
	iTemp = SingleToBuffer(dTemp*100.0, 14);
	SaveEMEM(01572, iTemp);

	//TAZ - Azimuth of landmark 1 at launch site
	dTemp = 309.6 - 360.0; //Degrees
	iTemp = SingleToBuffer(dTemp / 180.0, 0);
	SaveEMEM(01346, iTemp);

	//TAZ+1 - Azimuth of landmark 2 at launch site
	dTemp = 270.0 - 360.0; //Degrees
	iTemp = SingleToBuffer(dTemp / 180.0, 0);
	SaveEMEM(01347, iTemp);

	//TEL - Elevation of landmark 1 at launch site
	dTemp = -0.1; //Degrees
	iTemp = SingleToBuffer(dTemp / 180.0, 0);
	SaveEMEM(01350, iTemp);

	//TEL+1 - Elevation of landmark 2 at launch site
	dTemp = 0.0; //Degrees
	iTemp = SingleToBuffer(dTemp / 180.0, 0);
	SaveEMEM(01351, iTemp);

	//TATLAN1 - Nominal flight time to Atlantic target
	dTemp = BLOCKI.T_ATL; //Seconds
	DoubleToBuffer(dTemp*100.0, 28, iTemp, iTemp2);
	SaveEMEM(01617, iTemp);
	SaveEMEM(01620, iTemp2);

	//TPACIF1 - Nominal flight time to pacific target
	dTemp = BLOCKI.T_PAC; //Seconds
	DoubleToBuffer(dTemp*100.0, 28, iTemp, iTemp2);
	SaveEMEM(01627, iTemp);
	SaveEMEM(01630, iTemp2);

	//MISSION CONTROL PROGRAM

	//TDECAY - effective thrust decay time
	dTemp = -0.49; //seconds
	iTemp = SingleToBuffer(dTemp * 100, 14);
	SaveEMEM(01560, iTemp);

	//DELTAT - Value of computing interval
	dTemp = 2.0; //Seconds
	DoubleToBuffer(dTemp*100.0, 9, iTemp, iTemp2);
	SaveEMEM(01027, iTemp);
	SaveEMEM(01030, iTemp2);

	//ESQ(VR) - eccentricity squared for SPS1 burn
	DoubleToBuffer(BLOCKI.e_SPS1*BLOCKI.e_SPS1, 4, iTemp, iTemp2);
	SaveEMEM(01546, iTemp);
	SaveEMEM(01547, iTemp2);

	//ESQ(VR)+2 - eccentricity squared for SPS2 burn
	DoubleToBuffer(BLOCKI.e_SPS2*BLOCKI.e_SPS2, 4, iTemp, iTemp2);
	SaveEMEM(01550, iTemp);
	SaveEMEM(01551, iTemp2);

	//SEMILAT - semi-latus rectum for SPS1 burn
	dTemp = BLOCKI.a_SPS1*(1.0 - pow(BLOCKI.e_SPS1, 2));
	DoubleToBuffer(dTemp, 27, iTemp, iTemp2);
	SaveEMEM(01552, iTemp);
	SaveEMEM(01553, iTemp2);

	//SEMILAT+2 - semi-latus rectum for SPS2 burn
	dTemp = BLOCKI.a_SPS2*(1.0 - pow(BLOCKI.e_SPS2, 2));
	DoubleToBuffer(dTemp, 27, iTemp, iTemp2);
	SaveEMEM(01554, iTemp);
	SaveEMEM(01555, iTemp2);
}

void AGCPadloadGenerator::CoronaDefaults()
{
	BlockIDefaults();

	//POLYCOEF
	int tempaddr = 01575;
	for (int N = 0; N < 7; N++)
	{
		dTemp = BLOCKI.POLYCOFF[N] / (180.0*pow(100.0, N))*pow(2, 14 * N - 4);
		DoubleToBuffer(dTemp, 0, iTemp, iTemp2);
		SaveEMEM(tempaddr, iTemp);
		SaveEMEM(tempaddr + 1, iTemp2);
		tempaddr += 2;
	}

	//NSHIFT - Average G routine scaling constant
	SaveEMEM(01040, 077772);

	//XSHIFT - Average G routine scaling constant
	SaveEMEM(01041, 011);

	//CALCG - Preload CADR of CALCGEAR
	SaveEMEM(01042, 061652);

	MATRIX3 REFS = SolariumEarthFixedToSM(PadLat*RAD, PadLong*RAD, LaunchAzimuth*RAD);

	//RTATLAN1 - Post-LET abort target vector at lif-off + TATLAN1 in IMU coordinates assuming the platform goes inertial
	VECTOR3 RATL = unit(r_from_latlong(BLOCKI.lat_ATL*RAD, BLOCKI.lng_ATL*RAD + w_Earth * BLOCKI.T_ATL, 0.0, 0, 0));
	VECTOR3 RN = mul(REFS, RATL);
	DoubleToBuffer(RN.x, 1, iTemp, iTemp2);
	SaveEMEM(01621, iTemp);
	SaveEMEM(01622, iTemp2);
	DoubleToBuffer(RN.y, 1, iTemp, iTemp2);
	SaveEMEM(01623, iTemp);
	SaveEMEM(01624, iTemp2);
	DoubleToBuffer(RN.z, 1, iTemp, iTemp2);
	SaveEMEM(01625, iTemp);
	SaveEMEM(01626, iTemp2);

	//RTPACIF1 - Post-LET abort target vector at lif-off + TPACIF1 in IMU coordinates assuming the platform goes inertial
	VECTOR3 RPAC = unit(r_from_latlong(BLOCKI.lat_PAC*RAD, BLOCKI.lng_PAC*RAD + w_Earth * BLOCKI.T_PAC, 0.0, 0, 0));
	RN = mul(REFS, RPAC);
	DoubleToBuffer(RN.x, 1, iTemp, iTemp2);
	SaveEMEM(01631, iTemp);
	SaveEMEM(01632, iTemp2);
	DoubleToBuffer(RN.y, 1, iTemp, iTemp2);
	SaveEMEM(01633, iTemp);
	SaveEMEM(01634, iTemp2);
	DoubleToBuffer(RN.z, 1, iTemp, iTemp2);
	SaveEMEM(01635, iTemp);
	SaveEMEM(01636, iTemp2);

	//UNITW

	VECTOR3 UNITW = mul(REFS, _V(0, 0, 1));
	DoubleToBuffer(UNITW.x, 1, iTemp, iTemp2);
	SaveEMEM(01043, iTemp);
	SaveEMEM(01044, iTemp2);
	DoubleToBuffer(UNITW.y, 1, iTemp, iTemp2);
	SaveEMEM(01045, iTemp);
	SaveEMEM(01046, iTemp2);
	DoubleToBuffer(UNITW.z, 1, iTemp, iTemp2);
	SaveEMEM(01047, iTemp);
	SaveEMEM(01050, iTemp2);

	//RN - Position vector at GRR
	VECTOR3 R_L = r_from_latlong(PadLat*RAD, PadLong*RAD, PadAlt, BODY_EARTH, 0);

	RN = mul(REFS, R_L);
	DoubleToBuffer(RN.x, 24, iTemp, iTemp2);
	SaveEMEM(0765, iTemp);
	SaveEMEM(0766, iTemp2);
	DoubleToBuffer(RN.y, 24, iTemp, iTemp2);
	SaveEMEM(0767, iTemp);
	SaveEMEM(0770, iTemp2);
	DoubleToBuffer(RN.z, 24, iTemp, iTemp2);
	SaveEMEM(0771, iTemp);
	SaveEMEM(0772, iTemp2);

	//TCOAST - Time from SPS1 cutoff to maneuver to SPS2 attitude
	dTemp = 3163.67; //Seconds
	DoubleToBuffer(dTemp*100.0, 28, iTemp, iTemp2);
	SaveEMEM(01556, iTemp);
	SaveEMEM(01557, iTemp2);
}

void AGCPadloadGenerator::SolariumDefaults()
{
	BlockIDefaults();

	//POLYCOEF
	int tempaddr = 01575;
	for (int N = 0; N < 7; N++)
	{
		dTemp = BLOCKI.POLYCOFF[N] / (360.0*pow(100.0, N))*pow(2, 14 * N - 4);
		DoubleToBuffer(dTemp, 0, iTemp, iTemp2);
		SaveEMEM(tempaddr, iTemp);
		SaveEMEM(tempaddr + 1, iTemp2);
		tempaddr += 2;
	}

	//NSHIFT - Average G routine scaling constant
	SaveEMEM(01040, 077772);

	//XSHIFT - Average G routine scaling constant
	SaveEMEM(01041, 011);

	//MAXROLL - FInal roll angle minus initial roll angle
	dTemp = 18.0; //Degrees
	DoubleToBuffer(dTemp / 360.0, 0, iTemp, iTemp2);
	SaveEMEM(01702, iTemp);
	SaveEMEM(01703, iTemp2);

	//1/RLLRTE - One over desired roll rate
	dTemp = 1.0; //Degrees per second
	DoubleToBuffer(dTemp * 100.0 *360.0, 28, iTemp, iTemp2);
	SaveEMEM(01700, iTemp);
	SaveEMEM(01701, iTemp2);

	MATRIX3 REFS = SolariumEarthFixedToSM(PadLat*RAD, PadLong*RAD, LaunchAzimuth*RAD);

	//RTATLAN1 - Post-LET abort target vector at lif-off + TATLAN1 in IMU coordinates assuming the platform goes inertial
	VECTOR3 RATL = unit(r_from_latlong(BLOCKI.lat_ATL*RAD, BLOCKI.lng_ATL*RAD + w_Earth * BLOCKI.T_ATL, 0.0, 0, 0));
	VECTOR3 RN = mul(REFS, RATL);
	DoubleToBuffer(RN.x, 1, iTemp, iTemp2);
	SaveEMEM(01621, iTemp);
	SaveEMEM(01622, iTemp2);
	DoubleToBuffer(RN.y, 1, iTemp, iTemp2);
	SaveEMEM(01623, iTemp);
	SaveEMEM(01624, iTemp2);
	DoubleToBuffer(RN.z, 1, iTemp, iTemp2);
	SaveEMEM(01625, iTemp);
	SaveEMEM(01626, iTemp2);

	//RTPACIF1 - Post-LET abort target vector at lif-off + TPACIF1 in IMU coordinates assuming the platform goes inertial
	VECTOR3 RPAC = unit(r_from_latlong(BLOCKI.lat_PAC*RAD, BLOCKI.lng_PAC*RAD + w_Earth * BLOCKI.T_PAC, 0.0, 0, 0));
	RN = mul(REFS, RPAC);
	DoubleToBuffer(RN.x, 1, iTemp, iTemp2);
	SaveEMEM(01631, iTemp);
	SaveEMEM(01632, iTemp2);
	DoubleToBuffer(RN.y, 1, iTemp, iTemp2);
	SaveEMEM(01633, iTemp);
	SaveEMEM(01634, iTemp2);
	DoubleToBuffer(RN.z, 1, iTemp, iTemp2);
	SaveEMEM(01635, iTemp);
	SaveEMEM(01636, iTemp2);

	//UNITW
	
	VECTOR3 UNITW = mul(REFS, _V(0, 0, 1));
	DoubleToBuffer(UNITW.x, 1, iTemp, iTemp2);
	SaveEMEM(01043, iTemp);
	SaveEMEM(01044, iTemp2);
	DoubleToBuffer(UNITW.y, 1, iTemp, iTemp2);
	SaveEMEM(01045, iTemp);
	SaveEMEM(01046, iTemp2);
	DoubleToBuffer(UNITW.z, 1, iTemp, iTemp2);
	SaveEMEM(01047, iTemp);
	SaveEMEM(01050, iTemp2);

	//RN - Position vector at GRR
	VECTOR3 R_L = r_from_latlong(PadLat*RAD, PadLong*RAD, PadAlt, BODY_EARTH, 0);
	
	RN = mul(REFS, R_L);
	DoubleToBuffer(RN.x, 25, iTemp, iTemp2);
	SaveEMEM(0765, iTemp);
	SaveEMEM(0766, iTemp2);
	DoubleToBuffer(RN.y, 25, iTemp, iTemp2);
	SaveEMEM(0767, iTemp);
	SaveEMEM(0770, iTemp2);
	DoubleToBuffer(RN.z, 25, iTemp, iTemp2);
	SaveEMEM(0771, iTemp);
	SaveEMEM(0772, iTemp2);

	//MISSION CONTROL PROGRAM	

	//TFFMIN - Time-to-free-fall at which to ??? SPS2 ignition in 2 minutes
	dTemp = 687.0; //Seconds
	DoubleToBuffer(dTemp*100.0, 28, iTemp, iTemp2);
	SaveEMEM(01676, iTemp);
	SaveEMEM(01677, iTemp2);

	//TFFNOM - Value of TFF to use to compute time-off coast if the TT is not computable
	dTemp = 16200.0; //Seconds
	DoubleToBuffer(dTemp*100.0, 28, iTemp, iTemp2);
	SaveEMEM(01720, iTemp);
	SaveEMEM(01721, iTemp2);

	//CGY - SPS1 C.G. rotation around S/C Y-axis
	dTemp = 0.0229; //radians
	DoubleToBuffer(dTemp, 0, iTemp, iTemp2);
	SaveEMEM(01704, iTemp);
	SaveEMEM(01705, iTemp2);

	//CGY+2 - SPS2 C.G. rotation around S/C Y-axis
	dTemp = 0.0297; //radians
	DoubleToBuffer(dTemp, 0, iTemp, iTemp2);
	SaveEMEM(01706, iTemp);
	SaveEMEM(01707, iTemp2);

	//CGZ - SPS1 C.G. rotation around S/C Z-axis
	dTemp = 0.0833; //radians
	DoubleToBuffer(dTemp, 0, iTemp, iTemp2);
	SaveEMEM(01710, iTemp);
	SaveEMEM(01711, iTemp2);

	//CGZ+2 - SPS2 C.G. rotation around S/C Z-axis
	dTemp = 0.0995; //radians
	DoubleToBuffer(dTemp, 0, iTemp, iTemp2);
	SaveEMEM(01712, iTemp);
	SaveEMEM(01713, iTemp2);

	//ATDT - SPS1 integratied initial thrust acceleration magnitude
	dTemp = 9.589; // m/s
	DoubleToBuffer(dTemp / 100.0, 5, iTemp, iTemp2);
	SaveEMEM(01714, iTemp);
	SaveEMEM(01715, iTemp2);

	//ATDT+2 - SPS2 integratied initial thrust acceleration magnitude
	dTemp = 30.48; // m/s
	DoubleToBuffer(dTemp / 100.0, 5, iTemp, iTemp2);
	SaveEMEM(01716, iTemp);
	SaveEMEM(01717, iTemp2);

	//S2SWITCH - switch to re-compute SPS2 burn attitude
	SaveEMEM(01722, 0);

	//REFSWITCH - switch to force 280k ft free fall reference
	SaveEMEM(01723, 0);

	//REDOSPS1 - switch to repeat SPS1 at SPS2 ignition
	SaveEMEM(01724, 0);

	//ANGLEX, Y, Z - desired cold soak gimbal angles
	dTemp = 53.23; //degrees
	iTemp = SingleToBuffer(dTemp / 180.0, 0);
	SaveEMEM(01673, iTemp);

	dTemp = 278.99 - 360.0; //degrees
	iTemp = SingleToBuffer(dTemp / 180.0, 0);
	SaveEMEM(01674, iTemp);

	dTemp = 13.19; //degrees
	iTemp = SingleToBuffer(dTemp / 180.0, 0);
	SaveEMEM(01675, iTemp);

	//UPTIME - Time to incorporate 1st RVT update
	SaveEMEM(01671, 037777);
	SaveEMEM(01672, 037777);
}

void AGCPadloadGenerator::Skylark048Padload()
{
	//FLAGWRD1
	SaveEMEM(075, 0);
	//FLAGWRD3
	SaveEMEM(077, 040000);
	//FLAGWRD10
	SaveEMEM(0106, 0);
	//C31FLWRD
	SaveEMEM(0373, 0);
	//N26/PRI
	dTemp = 0.0; //minutes
	iTemp = SingleToBuffer(dTemp*60.0*100.0, 14);
	SaveEMEM(01016, iTemp);
	//N26/2CAD
	SaveEMEM(01017, 0);
	SaveEMEM(01020, 0);
	//PIPTIME
	double PIPTIME = (LaunchMJD - PrelaunchMJD)*8.64e6;
	DoubleToBuffer(PIPTIME, 28, iTemp, iTemp2);
	SaveEMEM(01035, iTemp);
	SaveEMEM(01036, iTemp2);
	//PGNCSALT
	dTemp = PadAlt;
	DoubleToBuffer(dTemp, 29, iTemp, iTemp2);
	SaveEMEM(01122, iTemp);
	SaveEMEM(01123, iTemp2);
	//PADLONG
	dTemp = PadLong / 360.0;
	DoubleToBuffer(dTemp, 0, iTemp, iTemp2);
	SaveEMEM(01124, iTemp);
	SaveEMEM(01125, iTemp2);
	//FIXTIME
	dTemp = 0.0; //sec
	DoubleToBuffer(dTemp*100.0, 28, iTemp, iTemp2);
	SaveEMEM(01333, iTemp);
	SaveEMEM(01334, iTemp2);
	//CDUCHKWD
	iTemp = SingleToBuffer(CDUCHKWD * 100.0, 14);
	SaveEMEM(01356, iTemp);

	IMUCompensation(true, false);

	//EMDOT
	dTemp = 64.89; //lbs/second
	iTemp = SingleToBuffer(dTemp*LBS2KG / 100.0, 3);
	SaveEMEM(01741, iTemp);

	//RVAR
	dTemp = 0.0; //PERC ERR
	DoubleToBuffer(dTemp, -16, iTemp, iTemp2);
	SaveEMEM(01742, iTemp);
	SaveEMEM(01743, iTemp2);

	//RVARMIN
	dTemp = 40000.0; //ft^2
	TripleToBuffer(-dTemp * pow(FT2M, 2), 40, iTemp, iTemp2, iTemp3);
	SaveEMEM(01744, iTemp);
	SaveEMEM(01745, iTemp2);
	SaveEMEM(01746, iTemp3);

	//EIMP1SEC
	dTemp = 19036.0; // lbf
	iTemp = SingleToBuffer(dTemp*LBF2N / 100.0, 14);
	SaveEMEM(01747, iTemp);

	//EFIMP01
	dTemp = 24028.0; // lbf
	iTemp = SingleToBuffer(dTemp*LBF2N / 100.0, 14);
	SaveEMEM(01750, iTemp);

	//EFIMP16
	dTemp = 20324.0; // lbf
	iTemp = SingleToBuffer(dTemp*LBF2N / 100.0, 14);
	SaveEMEM(01751, iTemp);

	//LADPAD
	dTemp = BLOCKII.LADPAD;
	iTemp = SingleToBuffer(dTemp, 0);
	SaveEMEM(01752, iTemp);
	//LODPAD
	dTemp = BLOCKII.LODPAD;
	iTemp = SingleToBuffer(dTemp, 0);
	SaveEMEM(01753, iTemp);

	//DVTHRESH
	iTemp = SingleToBuffer(DVTHRESH / 100.0, -2);
	SaveEMEM(01754, iTemp);

	//ALTVAR
	iTemp = SingleToBuffer(ALTVAR, -16);
	SaveEMEM(01755, iTemp);

	//WRENDPOS
	iTemp = SingleToBuffer(BLOCKII.WRENDPOS*0.3048, 19);
	SaveEMEM(02000, iTemp);
	//WRENDVEL
	iTemp = SingleToBuffer(BLOCKII.WRENDVEL*0.3048 / 100.0, 0);
	SaveEMEM(02001, iTemp);
	//RMAX
	iTemp = SingleToBuffer(BLOCKII.RMAX*0.3048, 19);
	SaveEMEM(02002, iTemp);
	//VMAX
	iTemp = SingleToBuffer(BLOCKII.VMAX*0.3048 / 100.0, 7);
	SaveEMEM(02003, iTemp);

	//EMSALT
	DoubleToBuffer(EMSALT*0.3048, 29, iTemp, iTemp2);
	SaveEMEM(02004, iTemp);
	SaveEMEM(02005, iTemp2);

	//TCS
	dTemp = 2220.0; //sec
	DoubleToBuffer(dTemp*100.0, 28, iTemp, iTemp2);
	SaveEMEM(02040, iTemp);
	SaveEMEM(02041, iTemp2);

	//EPS1
	dTemp = 0.0025; //deg
	DoubleToBuffer(dTemp / 360.0, 0, iTemp, iTemp2);
	SaveEMEM(02042, iTemp);
	SaveEMEM(02043, iTemp2);

	//DHNCC
	dTemp = 20.0; //NM
	DoubleToBuffer(dTemp*1852.0, 29, iTemp, iTemp2);
	SaveEMEM(02044, iTemp);
	SaveEMEM(02045, iTemp2);

	//EPS2
	dTemp = 1000.0; //ft
	DoubleToBuffer(dTemp*FT2M, 29, iTemp, iTemp2);
	SaveEMEM(02046, iTemp);
	SaveEMEM(02047, iTemp2);

	//DELH1
	dTemp = 10.0; //NM
	DoubleToBuffer(dTemp*1852.0, 29, iTemp, iTemp2);
	SaveEMEM(02050, iTemp);
	SaveEMEM(02051, iTemp2);

	//NH1
	iTemp = SingleToBuffer(0.5, 5);
	SaveEMEM(02052, iTemp);

	//WRDTIME
	dTemp = 2400.0; //sec
	iTemp = SingleToBuffer(dTemp*100.0, 28);
	SaveEMEM(02240, iTemp);

	//MINBLKTM
	dTemp = 328.8; //sec
	iTemp = SingleToBuffer(dTemp*100.0, 28);
	SaveEMEM(02241, iTemp);

	//TBEFCCMP
	dTemp = 822.0; //sec
	iTemp = SingleToBuffer(dTemp*100.0, 28);
	SaveEMEM(02242, iTemp);

	//BRNBLKTM
	dTemp = 822.0; //sec
	iTemp = SingleToBuffer(dTemp*100.0, 28);
	SaveEMEM(02243, iTemp);

	//MAXWTIME
	dTemp = 3616.8; //sec
	iTemp = SingleToBuffer(dTemp*100.0, 28);
	SaveEMEM(02244, iTemp);

	//FINCMPTM
	dTemp = 493.2; //sec
	iTemp = SingleToBuffer(dTemp*100.0, 28);
	SaveEMEM(02245, iTemp);

	//INTVAR
	dTemp = 3600.0; //m^2
	iTemp = SingleToBuffer(dTemp, 15);
	SaveEMEM(02246, iTemp);

	//NBOA(YZ)
	dTemp = 0.0;
	DoubleToBuffer(dTemp, 1, iTemp, iTemp2);
	SaveEMEM(02255, iTemp);
	SaveEMEM(02256, iTemp2);

	//NBOA(YZ)+2
	dTemp = cos(35.0*RAD);
	DoubleToBuffer(dTemp, 1, iTemp, iTemp2);
	SaveEMEM(02257, iTemp);
	SaveEMEM(02260, iTemp2);

	//NBOA(YZ)+4
	dTemp = -sin(35.0*RAD);
	DoubleToBuffer(dTemp, 1, iTemp, iTemp2);
	SaveEMEM(02261, iTemp);
	SaveEMEM(02262, iTemp2);

	//NBOA(YZ)+6
	dTemp = 0.0;
	DoubleToBuffer(dTemp, 1, iTemp, iTemp2);
	SaveEMEM(02263, iTemp);
	SaveEMEM(02264, iTemp2);

	//NBOA(YZ)+8
	dTemp = -sin(35.0*RAD);
	DoubleToBuffer(dTemp, 1, iTemp, iTemp2);
	SaveEMEM(02265, iTemp);
	SaveEMEM(02266, iTemp2);

	//NBOA(YZ)+10
	dTemp = -cos(35.0*RAD);
	DoubleToBuffer(dTemp, 1, iTemp, iTemp2);
	SaveEMEM(02267, iTemp);
	SaveEMEM(02270, iTemp2);

	//RTRIDOT
	dTemp = 7.32e-5; //meters/sec^3
	DoubleToBuffer(dTemp/pow(100.0, 3), -31, iTemp, iTemp2);
	SaveEMEM(02334, iTemp);
	SaveEMEM(02335, iTemp2);

	//1/2ALPHA
	dTemp = CMCDATA.C12ALPHA*100.0*360.0;
	DoubleToBuffer(dTemp, 19, iTemp, iTemp2);
	SaveEMEM(02373, iTemp);
	SaveEMEM(02374, iTemp2);

	//AZIMUTH
	dTemp = PadAzimuth / 360.0;
	DoubleToBuffer(dTemp, 0, iTemp, iTemp2);
	SaveEMEM(02400, iTemp);
	SaveEMEM(02401, iTemp2);

	//LATITUDE
	dTemp = PadLat / 360.0;
	DoubleToBuffer(dTemp, 0, iTemp, iTemp2);
	SaveEMEM(02402, iTemp);
	SaveEMEM(02403, iTemp2);

	//TAZEL
	dTemp = TAZEL[0] / 360.0;
	iTemp = SingleToBuffer(dTemp, -1);
	if (iTemp > 037777) iTemp++; //Two's complement
	SaveEMEM(02432, iTemp);
	dTemp = TAZEL[1] / 360.0;
	iTemp = SingleToBuffer(dTemp, -2);
	SaveEMEM(02433, iTemp);
	dTemp = TAZEL[2] / 360.0;
	iTemp = SingleToBuffer(dTemp, -1);
	if (iTemp > 037777) iTemp++; //Two's complement
	SaveEMEM(02434, iTemp);
	dTemp = TAZEL[3] / 360.0;
	iTemp = SingleToBuffer(dTemp, -2);
	SaveEMEM(02435, iTemp);

	//LAUNCHAZ
	dTemp = LaunchAzimuth / 360.0;
	DoubleToBuffer(dTemp, 0, iTemp, iTemp2);
	SaveEMEM(02633, iTemp);
	SaveEMEM(02634, iTemp2);

	//ETDECAY
	dTemp = 0.6*100.0;
	iTemp = SingleToBuffer(dTemp, 14);
	SaveEMEM(03000, iTemp);

	//EKPRIME
	SaveEMEM(03001, 0123);
	SaveEMEM(03002, 0175);
	//EKTLX
	SaveEMEM(03003, 017433);
	SaveEMEM(03004, 04500);
	SaveEMEM(03005, 0334);
	//EREPFRAC
	SaveEMEM(03006, 01000);
	SaveEMEM(03007, 0232);
	//PACTOFF
	iTemp = SingleToBuffer(BLOCKII.PACTOFF / 0.023725, 14);
	SaveEMEM(03010, iTemp);
	//YACTOFF
	iTemp = SingleToBuffer(BLOCKII.YACTOFF / 0.023725, 14);
	SaveEMEM(03011, iTemp);
	//HBN10
	SaveEMEM(03012, 037777);
	//HBN11/2
	SaveEMEM(03013, 0);
	//HBN12
	SaveEMEM(03014, 0);
	//HBD11/2
	SaveEMEM(03015, 054360);
	//HBD12
	SaveEMEM(03016, 021075);
	//HBN20
	SaveEMEM(03017, 037777);
	//HBN21/2
	SaveEMEM(03020, 060465);
	//HBN22
	SaveEMEM(03021, 0);
	//HBD21/2
	SaveEMEM(03022, 054360);
	//HBD22
	SaveEMEM(03023, 021075);
	//HBN30
	SaveEMEM(03024, 037777);
	//HBN31/2
	SaveEMEM(03025, 057142);
	//HBN32
	SaveEMEM(03026, 033106);
	//HBD31/2
	SaveEMEM(03027, 050741);
	//HBD32
	SaveEMEM(03030, 031162);
	//SLOPE2+0
	SaveEMEM(03031, 012173);
	//SLOPE2+1
	SaveEMEM(03032, 07534);
	//WLH/SLOP+0
	SaveEMEM(03033, 0114);
	//WLH/SLOP+1
	SaveEMEM(03034, 057);
	//WL-H/SLP+0
	SaveEMEM(03035, 056);
	//WL-H/SLP+1
	SaveEMEM(03036, 017);

	//ECP
	iTemp = SingleToBuffer(CMCDATA.ECP, 0);
	SaveEMEM(03037, iTemp);

	//ECYW
	iTemp = SingleToBuffer(CMCDATA.ECYW, 0);
	SaveEMEM(03040, iTemp);

	//ALPHAP
	iTemp = SingleToBuffer(CMCDATA.ALPHAP, 0);
	SaveEMEM(03041, iTemp);

	//ALPHAYW
	iTemp = SingleToBuffer(CMCDATA.ALPHAYW, 0);
	SaveEMEM(03042, iTemp);

	//KMJDCKD
	iTemp = SingleToBuffer(CMCDATA.KMJDCKD / 1000.0 / 360.0, -13);
	SaveEMEM(03043, iTemp);

	//KMJ1DCKD
	iTemp = SingleToBuffer(CMCDATA.KMJ1DCKD / 1000.0 / 360.0, -13);
	SaveEMEM(03044, iTemp);
	
	//KMJ2DCKD
	iTemp = SingleToBuffer(CMCDATA.KMJ2DCKD / 1000.0 / 360.0, -13);
	SaveEMEM(03045, iTemp);

	//J/MDCKD
	iTemp = SingleToBuffer(CMCDATA.JMDCKD*1000.0*360.0, 27);
	SaveEMEM(03046, iTemp);

	//J/M1DCKD
	iTemp = SingleToBuffer(CMCDATA.JM1DCKD*1000.0*360.0, 27);
	SaveEMEM(03047, iTemp);

	//J/M2DCKD
	iTemp = SingleToBuffer(CMCDATA.JM2DCKD*1000.0*360.0, 27);
	SaveEMEM(03050, iTemp);

	//CLDKDELT
	SaveEMEM(03055, 0247);
	//DAPDATR3
	SaveEMEM(03070, 011111);
	//CH5FAIL
	SaveEMEM(03071, 0146);
	//CH6FAIL
	SaveEMEM(03072, CMCDATA.CH6FAIL);
	//DKRATE
	iTemp = SingleToBuffer(CMCDATA.DKRATE / 10.0 / 360.0, -9);
	SaveEMEM(03073, iTemp);
	//DKDB
	SaveEMEM(03074, 056);
	//WHICHDAP
	SaveEMEM(03075, 0);
	//WHICHX2
	SaveEMEM(03076, 0);
	//DAPDATR1
	SaveEMEM(03114, 031102);
	//DAPDATR2
	SaveEMEM(03115, 01111);

	//LEMMASS
	iTemp = SingleToBuffer(BLOCKII.LMMass*LBS2KG, 16);
	SaveEMEM(03121, iTemp);
	//CSMMASS
	iTemp = SingleToBuffer(BLOCKII.CSMMass*LBS2KG, 16);
	SaveEMEM(03122, iTemp);

	//POLYNUM
	SavePOLYNUM(03312);

	//SATRLRT
	SaveEMEM(03331, 0);
	SaveEMEM(03332, 016441);
	//RPSTART
	dTemp = BLOCKII.RPSTART*100.0;
	iTemp = SingleToBuffer(dTemp, 14);
	SaveEMEM(03333, iTemp);
	//POLYSTOP
	dTemp = -BLOCKII.POLYSTOP*100.0;
	iTemp = SingleToBuffer(dTemp, 14);
	SaveEMEM(03334, iTemp);
	//SATRATE
	SaveEMEM(03341, 0);
	SaveEMEM(03342, 0344);
	SaveEMEM(03343, 077433);
	SaveEMEM(03344, 0);
	//SATSCALE
	SaveEMEM(03351, 010000); //0.3 V/DEG

	//ALFAPAD
	dTemp = BLOCKII.ALFAPAD / 360.0;
	iTemp = SingleToBuffer(dTemp, -1);
	SaveEMEM(03377, iTemp);

	//LAT(SPL)
	dTemp = BLOCKII.LAT_SPL / 360.0;
	DoubleToBuffer(dTemp, 0, iTemp, iTemp2);
	SaveEMEM(03400, iTemp);
	SaveEMEM(03401, iTemp2);
	//LNG(SPL)
	dTemp = BLOCKII.LNG_SPL / 360.0;
	DoubleToBuffer(dTemp, 0, iTemp, iTemp2);
	SaveEMEM(03402, iTemp);
	SaveEMEM(03403, iTemp2);

	//SBDELT
	dTemp = 0.81*100.0;
	iTemp = SingleToBuffer(dTemp, 14);
	SaveEMEM(03424, iTemp);

	//ELEV
	dTemp = 27.0; //deg
	DoubleToBuffer(dTemp / 360.0, 0, iTemp, iTemp2);
	SaveEMEM(03640, iTemp);
	SaveEMEM(03641, iTemp2);

	//TTPI
	dTemp = CMCDATA.TTPI*100.0;
	DoubleToBuffer(dTemp, 28, iTemp, iTemp2);
	SaveEMEM(03642, iTemp);
	SaveEMEM(03643, iTemp2);
}

void AGCPadloadGenerator::SkylarkSolarEphemeris(double TC, double T0)
{
	//t_c: julian ephemeris date at midpoint of the time interval over which the solar position approximation is desired
	//T0: julian ephemeris date midnight July 1 before launch

	double T, fact, NEM, e_EM, W_0EM, M_0EM, LOSO, CARG, CMOD;

	T = (TC - 2433282.5) / 36525.0;
	fact = 1.0 + 2.852e-8; //Converts ephemeris time rate to universal time rate
	NEM = fact * 0.9856002628;

	//Calculate mean eccentricity
	e_EM = 0.0167301085 - 4.1926e-5*T - 1.26e-7*T*T;
	//Calculate longitude of perihelion
	W_0EM = 102.08053 + (0.32328 / 36525.0)*(T0 - 2433282.5) + 1.5e-4*T*T;
	W_0EM = fmod(W_0EM, 360.0);
	W_0EM *= RAD;
	//Calculate mean anomaly
	M_0EM = 358.000682 + NEM * (T0 - 2433282.5) - 1.55e-4*pow(T, 2) - 3.3333e-6*pow(T, 3);
	M_0EM = fmod(M_0EM, 360.0);
	M_0EM *= RAD;

	LOSO = W_0EM + M_0EM - PI;
	CARG = M_0EM; // PHASEC
	CMOD = -(2.0 * e_EM - 0.25*pow(e_EM, 3)); //C

	//LOSO
	dTemp = LOSO;
	DoubleToBuffer(dTemp / PI2, 0, iTemp, iTemp2);
	SaveEMEM(02006, iTemp);
	SaveEMEM(02007, iTemp2);

	//CMOD
	dTemp = CMOD;
	DoubleToBuffer(dTemp / PI2, -1, iTemp, iTemp2);
	SaveEMEM(02010, iTemp);
	SaveEMEM(02011, iTemp2);

	//CARG
	dTemp = CARG;
	DoubleToBuffer(dTemp / PI2, 0, iTemp, iTemp2);
	SaveEMEM(02012, iTemp);
	SaveEMEM(02013, iTemp2);
}

void AGCPadloadGenerator::SkylarkCorrectionMatrix(double TC, double T0)
{
	//t_c: julian ephemeris date at midpoint of the time interval over which the correction matrix is desired
	//T0: julian ephemeris date midnight July 1 before launch

	MATRIX3 J2000, J2000_2, R, R3, R3_2, M, M_AGC;
	double w_E, MJD_C, MJD_0, A_Z, A_Z0;

	w_E = 7.29211514667e-5;   //Skylark

	MJD_C = JD2MJD(TC);
	MJD_0 = JD2MJD(T0);

	//J2000 = matrix converting from J2000 to BRCS
	J2000 = J2000EclToBRCS(1950);
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

	M = _M(cos(A_Z), sin(A_Z), 0.0, -sin(A_Z), cos(A_Z), 0.0, 0.0, 0.0, 1.0);
	M_AGC = tmat(mul(R3, M));
	
	//LMATRIX
	DoubleToBuffer(M_AGC.m11, 0, iTemp, iTemp2);
	SaveEMEM(02014, iTemp);
	SaveEMEM(02015, iTemp2);
	DoubleToBuffer(M_AGC.m12, 0, iTemp, iTemp2);
	SaveEMEM(02016, iTemp);
	SaveEMEM(02017, iTemp2);
	DoubleToBuffer(M_AGC.m13, 0, iTemp, iTemp2);
	SaveEMEM(02020, iTemp);
	SaveEMEM(02021, iTemp2);
	DoubleToBuffer(M_AGC.m21, 0, iTemp, iTemp2);
	SaveEMEM(02022, iTemp);
	SaveEMEM(02023, iTemp2);
	DoubleToBuffer(M_AGC.m22, 0, iTemp, iTemp2);
	SaveEMEM(02024, iTemp);
	SaveEMEM(02025, iTemp2);
	DoubleToBuffer(M_AGC.m23, 0, iTemp, iTemp2);
	SaveEMEM(02026, iTemp);
	SaveEMEM(02027, iTemp2);
	DoubleToBuffer(M_AGC.m31, 0, iTemp, iTemp2);
	SaveEMEM(02030, iTemp);
	SaveEMEM(02031, iTemp2);
	DoubleToBuffer(M_AGC.m32, 0, iTemp, iTemp2);
	SaveEMEM(02032, iTemp);
	SaveEMEM(02033, iTemp2);
	DoubleToBuffer(M_AGC.m33, 0, iTemp, iTemp2);
	SaveEMEM(02034, iTemp);
	SaveEMEM(02035, iTemp2);

	//AZ0
	dTemp = A_Z0;
	DoubleToBuffer(dTemp / PI2, 0, iTemp, iTemp2);
	SaveEMEM(02036, iTemp);
	SaveEMEM(02037, iTemp2);
}

MATRIX3 AGCPadloadGenerator::J2000EclToBRCSMJD(double mjd) const
{
	double t1 = (mjd - 51544.5) / 36525.0;
	double t2 = t1 * t1;
	double t3 = t2 * t1;

	t1 *= 4.848136811095359e-6;
	t2 *= 4.848136811095359e-6;
	t3 *= 4.848136811095359e-6;

	double i = 2004.3109*t1 - 0.42665*t2 - 0.041833*t3;
	double r = 2306.2181*t1 + 0.30188*t2 + 0.017998*t3;
	double L = 2306.2181*t1 + 1.09468*t2 + 0.018203*t3;

	double rot = -r - PI05;
	double lan = PI05 - L;
	double inc = i;
	double obl = 0.4090928023;

	return mul(mul(_MRz(rot), _MRx(inc)), mul(_MRz(lan), _MRx(-obl)));
}

MATRIX3 AGCPadloadGenerator::J2000EclToBRCS(int epoch) const
{
	//Calculate the rotation matrix between J2000 and mean Besselian of epoch coordinate systems 
	double MJD = MJDOfNBYEpoch(epoch);
	return J2000EclToBRCSMJD(MJD);
}

double AGCPadloadGenerator::TJUDAT(int Y, int M, int D) const
{
	//Year, Month, Day to Julian Date
	int Y_apo = Y - 1900;
	int TMM[] = { 0,31,59,90,120,151,181,212,243,273,304,334 };

	int Z = Y_apo / 4;
	if (Y_apo % 4 == 0)
	{
		Z = Z - 1;
		for (int i = 2; i < 12; i++)
		{
			TMM[i] += 1;
		}
	}
	return 2415020.5 + (double)(365 * Y_apo + Z + TMM[M - 1] + D - 1);
}

MATRIX3 AGCPadloadGenerator::GetRotationMatrix(int plan, double t) const
{
	double t0, T_p, L_0, e_rel, phi_0, T_s, e_ref, L_ref, L_rel, phi;
	MATRIX3 Rot1, Rot2, R_ref, Rot3, Rot4, R_rel, R_rot, R, Rot;

	if (plan == BODY_EARTH)
	{
		t0 = 51544.5;								//LAN_MJD, MJD of the LAN in the "beginning"
		T_p = -9413040.4;							//Precession Period
		L_0 = 0.00001553343;						//LAN in the "beginning"
		e_rel = 0.4090928023;						//Obliquity / axial tilt of the earth in radians
		phi_0 = 4.894942829;						//Sidereal Rotational Offset
		T_s = 86164.098904 / 24.0 / 60.0 / 60.0;	//Sidereal Rotational Period
		e_ref = 0;									//Precession Obliquity
		L_ref = 0;									//Precession LAN
	}
	else
	{
		t0 = 51544.5;							//LAN_MJD, MJD of the LAN in the "beginning"
		T_p = -6793.468728092782;				//Precession Period
		L_0 = 1.71817749;						//LAN in the "beginning"
		e_rel = 0.026699886264850;				//Obliquity / axial tilt of the earth in radians
		phi_0 = 4.769465382;					//Sidereal Rotational Offset
		T_s = 2360588.15 / 24.0 / 60.0 / 60.0;	//Sidereal Rotational Period
		e_ref = 7.259562816e-005;				//Precession Obliquity
		L_ref = 0.4643456618;					//Precession LAN
	}

	Rot1 = _M(cos(L_ref), 0, -sin(L_ref), 0, 1, 0, sin(L_ref), 0, cos(L_ref));
	Rot2 = _M(1, 0, 0, 0, cos(e_ref), -sin(e_ref), 0, sin(e_ref), cos(e_ref));
	R_ref = mul(Rot1, Rot2);
	L_rel = L_0 + PI2 * (t - t0) / T_p;
	Rot3 = _M(cos(L_rel), 0, -sin(L_rel), 0, 1, 0, sin(L_rel), 0, cos(L_rel));
	Rot4 = _M(1, 0, 0, 0, cos(e_rel), -sin(e_rel), 0, sin(e_rel), cos(e_rel));
	R_rel = mul(Rot3, Rot4);
	phi = phi_0 + PI2 * (t - t0) / T_s + (L_0 - L_rel)*cos(e_rel);
	R_rot = _M(cos(phi), 0, -sin(phi), 0, 1, 0, sin(phi), 0, cos(phi));
	Rot = mul(R_rel, R_rot);
	R = mul(R_ref, Rot);
	return R;
}

double AGCPadloadGenerator::MJDOfNBYEpoch(int epoch) const
{
	//Calculate MJD of Besselian epoch
	double C, DE, MJD, JD, T;
	int E, XN;
	const double A = 0.0929;
	const double B = 8640184.542;
	const double W1 = 1.720217954160054e-2;

	E = epoch;
	XN = (E - 1901) / 4;
	C = -86400.0*(double)(E - 1900) - 74.164;
	T = 2.0 * C / (-B - sqrt(B*B - 4.0 * A*C));
	DE = 36525.0*T - 365.0*(double)(E - 1900) + 0.5 - (double)XN;

	JD = TJUDAT(epoch, 1, 0);
	MJD = JD - 2400000.5 + DE;
	return MJD;
}

double AGCPadloadGenerator::MJD2JD(double MJD) const
{
	return MJD + 2400000.5;
}

double AGCPadloadGenerator::JD2MJD(double JD) const
{
	return JD - 2400000.5;
}