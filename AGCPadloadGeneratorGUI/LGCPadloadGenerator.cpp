// LGCPadloadGenerator.cpp
//

#include "pch.h"
#include "AGCPadloadGeneratorGUI.h"
#include "LGCPadloadGenerator.h"
#include "afxdialogex.h"

#define LGC_SUNDANCE306 0
#define LGC_LUMINARY069 1
#define LGC_LUMINARY069R2 2
#define LGC_LUMINARY099 3
#define LGC_LUMINARY116 4
#define LGC_LUMINARY131R1 5
#define LGC_LUMINARY178 6
#define LGC_LUMINARY210 7
#define LGC_ZERLINA056 8

// LGCPadloadGenerator-Dialog

IMPLEMENT_DYNAMIC(LGCPadloadGenerator, CDialogEx)

LGCPadloadGenerator::LGCPadloadGenerator(CWnd* pParent /*=nullptr*/)
	: CDialogEx(IDD_DIALOG2, pParent)
{

}

LGCPadloadGenerator::~LGCPadloadGenerator()
{
}

void LGCPadloadGenerator::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_COMBO1, RopeNameBox);
	DDX_Control(pDX, IDC_COMBO2, MissionBox);
	DDX_Control(pDX, IDC_EDIT1, LaunchMJDInput);
	DDX_Control(pDX, IDC_EDIT9, EphemerisSpanBox);
	DDX_Control(pDX, IDC_EDIT2, T504LMBox);
	DDX_Control(pDX, IDC_EDIT8, UNITWBox);
	DDX_Control(pDX, IDC_EDIT3, LSLatitudeBox);
	DDX_Control(pDX, IDC_EDIT4, LSLongitudeBox);
	DDX_Control(pDX, IDC_EDIT5, LSAltitudeBox);
	DDX_Control(pDX, IDC_EDIT6, LMMassBox);
	DDX_Control(pDX, IDC_EDIT7, CSMMassBox);
	DDX_Control(pDX, IDC_EDIT10, TotalMassBox);
	DDX_Control(pDX, IDC_CHECK1, DockedBox);
	DDX_Control(pDX, IDC_EDIT11, WRENDPOSBox);
	DDX_Control(pDX, IDC_EDIT12, WRENDVELBox);
	DDX_Control(pDX, IDC_EDIT13, WSHAFTBox);
	DDX_Control(pDX, IDC_EDIT14, WTRUNBox);
	DDX_Control(pDX, IDC_EDIT15, RMAXBox);
	DDX_Control(pDX, IDC_EDIT18, VMAXBox);
	DDX_Control(pDX, IDC_EDIT19, SHAFTVARBox);
	DDX_Control(pDX, IDC_EDIT20, TRUNVARBox);
	DDX_Control(pDX, IDC_EDIT21, WSURFPOSBox);
	DDX_Control(pDX, IDC_EDIT22, WSURFVELBox);
	DDX_Control(pDX, IDC_EDIT16, HIASCENTBox);
	DDX_Control(pDX, IDC_EDIT23, AGSKBox);
	DDX_Control(pDX, IDC_EDIT24, ROLLTIMEBox);
	DDX_Control(pDX, IDC_EDIT25, PITCHTIMEBox);
	DDX_Control(pDX, IDC_EDIT26, TLANDBox);
	DDX_Control(pDX, IDC_EDIT27, ABSC0Box);
	DDX_Control(pDX, IDC_EDIT28, ABSC1Box);
	DDX_Control(pDX, IDC_EDIT29, ABSC2Box);
	DDX_Control(pDX, IDC_EDIT30, ABSC3Box);
	DDX_Control(pDX, IDC_EDIT31, ABSC4Box);
	DDX_Control(pDX, IDC_EDIT32, SLOPE0Box);
	DDX_Control(pDX, IDC_EDIT33, SLOPE1Box);
	DDX_Control(pDX, IDC_EDIT34, SLOPE2Box);
	DDX_Control(pDX, IDC_EDIT35, SLOPE3Box);
	DDX_Control(pDX, IDC_EDIT36, SLOPE4Box);
	DDX_Control(pDX, IDC_EDIT37, IGNAOSQBox);
	DDX_Control(pDX, IDC_EDIT38, IGNAOSRBox);
	DDX_Control(pDX, IDC_EDIT39, VIGNBox);
	DDX_Control(pDX, IDC_EDIT40, RIGNXBox);
	DDX_Control(pDX, IDC_EDIT41, RIGNZBox);
	DDX_Control(pDX, IDC_EDIT42, KIGNXBox);
	DDX_Control(pDX, IDC_EDIT43, KIGNYBox);
	DDX_Control(pDX, IDC_EDIT44, KIGNVBox);
	DDX_Control(pDX, IDC_EDIT45, J1PARMBox);
	DDX_Control(pDX, IDC_EDIT46, K1PARMBox);
	DDX_Control(pDX, IDC_EDIT47, J2PARMBox);
	DDX_Control(pDX, IDC_EDIT48, K2PARMBox);
	DDX_Control(pDX, IDC_EDIT49, THETCRITBox);
	DDX_Control(pDX, IDC_EDIT50, RAMINBox);
	DDX_Control(pDX, IDC_EDIT51, DELTTFAPBox);
	DDX_Control(pDX, IDC_EDIT52, OutputBox);
	DDX_Control(pDX, IDC_EDIT56, PIOSDataSetBox);
	DDX_Control(pDX, IDC_CHECK2, R2ModelBox);
	DDX_Control(pDX, IDC_EDIT65, PBIASXBox);
	DDX_Control(pDX, IDC_EDIT66, PIPASCFXBox);
	DDX_Control(pDX, IDC_EDIT67, PBIASYBox);
	DDX_Control(pDX, IDC_EDIT68, PIPASCFYBox);
	DDX_Control(pDX, IDC_EDIT69, PBIASZBox);
	DDX_Control(pDX, IDC_EDIT70, PIPASCFZBox);
	DDX_Control(pDX, IDC_EDIT71, NBDXBox);
	DDX_Control(pDX, IDC_EDIT72, NBDYBox);
	DDX_Control(pDX, IDC_EDIT73, NBDZBox);
	DDX_Control(pDX, IDC_EDIT74, ADIAXBox);
	DDX_Control(pDX, IDC_EDIT75, ADIAYBox);
	DDX_Control(pDX, IDC_EDIT76, ADIAZBox);
	DDX_Control(pDX, IDC_EDIT77, ADSRAXBox);
	DDX_Control(pDX, IDC_EDIT78, ADSRAYBox);
	DDX_Control(pDX, IDC_EDIT79, ADSRAZBox);
}


BEGIN_MESSAGE_MAP(LGCPadloadGenerator, CDialogEx)
	ON_BN_CLICKED(IDOK, &LGCPadloadGenerator::OnBnClickedOk)
	ON_CBN_SELCHANGE(IDC_COMBO1, &LGCPadloadGenerator::OnCbnSelchangeCombo1)
	ON_CBN_SELCHANGE(IDC_COMBO2, &LGCPadloadGenerator::OnCbnSelchangeCombo2)
	ON_BN_CLICKED(IDCANCEL, &LGCPadloadGenerator::OnBnClickedCancel)
	ON_EN_CHANGE(IDC_EDIT6, &LGCPadloadGenerator::OnEnChangeEdit6)
	ON_EN_CHANGE(IDC_EDIT7, &LGCPadloadGenerator::OnEnChangeEdit7)
	ON_BN_CLICKED(IDC_CHECK1, &LGCPadloadGenerator::OnBnClickedCheck1)
END_MESSAGE_MAP()

BOOL LGCPadloadGenerator::PreTranslateMessage(MSG* pMsg)
{
	m_ToolTip.RelayEvent(pMsg);

	return CDialog::PreTranslateMessage(pMsg);
}

BOOL LGCPadloadGenerator::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	//Create the ToolTip control
	if (!m_ToolTip.Create(this))
	{
		TRACE0("Unable to create the ToolTip!");
	}
	else
	{
		m_ToolTip.AddTool(&PIOSDataSetBox, _T("Planetary Inertial Orientation Subroutine data from PIOSDataSets.txt"));

		m_ToolTip.AddTool(&UNITWBox, _T("Time in days from liftoff at which Earth rotations are exactly accurate"));
		m_ToolTip.AddTool(&T504LMBox, _T("Time in days from liftoff at which Moon rotations are exactly accurate"));

		m_ToolTip.AddTool(&WRENDPOSBox, _T("W matrix position initialization for rendezvous"));
		m_ToolTip.AddTool(&WRENDVELBox, _T("W matrix velocity initialization for rendezvous"));
		m_ToolTip.AddTool(&RMAXBox, _T("Maximum automatic rendezvous position update"));
		m_ToolTip.AddTool(&VMAXBox, _T("Maximum automatic rendezvous velocity update"));

		m_ToolTip.AddTool(&R2ModelBox, _T("R2 gravity model only supported by Open Orbiter"));
	}

	MissionBox.AddString(L"Manual");
	MissionBox.AddString(L"Apollo 9");
	MissionBox.AddString(L"Apollo 10");
	MissionBox.AddString(L"Apollo 11");
	MissionBox.AddString(L"Apollo 12");
	MissionBox.AddString(L"Apollo 13");
	MissionBox.AddString(L"Apollo 14");
	MissionBox.AddString(L"Apollo 15");
	MissionBox.AddString(L"Apollo 16");
	MissionBox.AddString(L"Apollo 17");
	MissionBox.SetCurSel(0);

	RopeNameBox.AddString(L"Sundance306");
	RopeNameBox.AddString(L"Luminary069");
	RopeNameBox.AddString(L"Luminary069R2");
	RopeNameBox.AddString(L"Luminary099");
	RopeNameBox.AddString(L"Luminary116");
	RopeNameBox.AddString(L"Luminary131R1");
	RopeNameBox.AddString(L"Luminary178");
	RopeNameBox.AddString(L"Luminary210");
	RopeNameBox.AddString(L"Zerlina56");
	RopeNameBox.SetCurSel(LGC_LUMINARY116);

	UpdateRopeSpecificEditFields();

	T504LMBox.SetWindowText(L"4.5");
	UNITWBox.SetWindowText(L"-0.5");
	LSLatitudeBox.SetWindowText(L"0.0");
	LSLongitudeBox.SetWindowText(L"0.0");
	LSAltitudeBox.SetWindowText(L"0.0");
	EphemerisSpanBox.SetWindowText(L"14.5");
	LaunchMJDInput.SetWindowText(L"0.0");
	LMMassBox.SetWindowText(L"33872.3");
	CSMMassBox.SetWindowText(L"37580.3");
	TotalMassBox.SetWindowText(L"33872.3");
	HIASCENTBox.SetWindowText(L"10900.0");

	WRENDPOSBox.SetWindowText(L"10000");
	WRENDVELBox.SetWindowText(L"10");
	WSHAFTBox.SetWindowText(L"15");
	WTRUNBox.SetWindowText(L"15");
	RMAXBox.SetWindowText(L"2000");
	VMAXBox.SetWindowText(L"2");
	WSURFPOSBox.SetWindowText(L"0");
	WSURFVELBox.SetWindowText(L"0");
	SHAFTVARBox.SetWindowText(L"1");
	TRUNVARBox.SetWindowText(L"1");

	AGSKBox.SetWindowText(L"100.0");

	ROLLTIMEBox.SetWindowText(L"6.0");
	PITCHTIMEBox.SetWindowText(L"6.0");

	TLANDBox.SetWindowText(L"0.0");

	ABSC0Box.SetWindowText(L"0.0");
	ABSC1Box.SetWindowText(L"0.0");
	ABSC2Box.SetWindowText(L"0.0");
	ABSC3Box.SetWindowText(L"0.0");
	ABSC4Box.SetWindowText(L"0.0");
	SLOPE0Box.SetWindowText(L"0.0");
	SLOPE1Box.SetWindowText(L"0.0");
	SLOPE2Box.SetWindowText(L"0.0");
	SLOPE3Box.SetWindowText(L"0.0");
	SLOPE4Box.SetWindowText(L"0.0");

	IGNAOSQBox.SetWindowText(L"0.0");
	IGNAOSRBox.SetWindowText(L"0.0");

	VIGNBox.SetWindowText(L"5546.447179");
	RIGNXBox.SetWindowText(L"-140345.7283");
	RIGNZBox.SetWindowText(L"-1464979.987");
	KIGNXBox.SetWindowText(L"-0.419");
	KIGNYBox.SetWindowText(L"-9.05e-7");
	KIGNVBox.SetWindowText(L"-470.0");

	J1PARMBox.SetWindowText(L"6032567.5");
	K1PARMBox.SetWindowText(L"-6.2726125e5");
	J2PARMBox.SetWindowText(L"6030470.0");
	K2PARMBox.SetWindowText(L"-3.1835146e5");
	THETCRITBox.SetWindowText(L"8.384852304");
	RAMINBox.SetWindowText(L"5.8768997e6");

	DELTTFAPBox.SetWindowText(L"-90.0");

	PBIASXBox.SetWindowText(L"0.0");
	PIPASCFXBox.SetWindowText(L"0.0");
	PBIASYBox.SetWindowText(L"0.0");
	PIPASCFYBox.SetWindowText(L"0.0");
	PBIASZBox.SetWindowText(L"0.0");
	PIPASCFZBox.SetWindowText(L"0.0");
	NBDXBox.SetWindowText(L"0.0");
	NBDYBox.SetWindowText(L"0.0");
	NBDZBox.SetWindowText(L"0.0");
	ADIAXBox.SetWindowText(L"0.0");
	ADIAYBox.SetWindowText(L"0.0");
	ADIAZBox.SetWindowText(L"0.0");
	ADSRAXBox.SetWindowText(L"0.0");
	ADSRAYBox.SetWindowText(L"0.0");
	ADSRAZBox.SetWindowText(L"0.0");

	return TRUE;
}

void LGCPadloadGenerator::OnBnClickedOk()
{
	CString string;

	agc.LaunchMJD = Utilities::Text2Double(&LaunchMJDInput);
	agc.T_504LM = Utilities::Text2Double(&T504LMBox);
	agc.T_UNITW = Utilities::Text2Double(&UNITWBox);
	agc.LSLat = Utilities::Text2Double(&LSLatitudeBox);
	agc.LSLng = Utilities::Text2Double(&LSLongitudeBox);
	agc.LSAlt = Utilities::Text2Double(&LSAltitudeBox);
	agc.EphemerisSpan = Utilities::Text2Double(&EphemerisSpanBox);
	agc.BLOCKII.LMMass = Utilities::Text2Double(&LMMassBox);
	agc.BLOCKII.CSMMass = Utilities::Text2Double(&CSMMassBox);
	agc.BLOCKII.TotalMass = Utilities::Text2Double(&TotalMassBox);
	agc.BLOCKII.HIASCENT = Utilities::Text2Double(&HIASCENTBox);
	agc.BLOCKII.WRENDPOS = Utilities::Text2Double(&WRENDPOSBox);
	agc.BLOCKII.WRENDVEL = Utilities::Text2Double(&WRENDVELBox);
	agc.BLOCKII.WSHAFT = Utilities::Text2Double(&WSHAFTBox);
	agc.BLOCKII.WTRUN = Utilities::Text2Double(&WTRUNBox);
	agc.BLOCKII.RMAX = Utilities::Text2Double(&RMAXBox);
	agc.BLOCKII.VMAX = Utilities::Text2Double(&VMAXBox);
	agc.BLOCKII.WSURFPOS = Utilities::Text2Double(&WSURFPOSBox);
	agc.BLOCKII.WSURFVEL = Utilities::Text2Double(&WSURFVELBox);
	agc.BLOCKII.SHAFTVAR = Utilities::Text2Double(&SHAFTVARBox);
	agc.BLOCKII.TRUNVAR = Utilities::Text2Double(&TRUNVARBox);
	agc.BLOCKII.AGSK = Utilities::Text2Double(&AGSKBox);
	agc.BLOCKII.ROLLTIME = Utilities::Text2Double(&ROLLTIMEBox);
	agc.BLOCKII.PITCHTIME = Utilities::Text2Double(&PITCHTIMEBox);
	agc.BLOCKII.TLAND = Utilities::Text2Double(&TLANDBox);
	agc.BLOCKII.ABSC[0] = Utilities::Text2Double(&ABSC0Box);
	agc.BLOCKII.ABSC[1] = Utilities::Text2Double(&ABSC1Box);
	agc.BLOCKII.ABSC[2] = Utilities::Text2Double(&ABSC2Box);
	agc.BLOCKII.ABSC[3] = Utilities::Text2Double(&ABSC3Box);
	agc.BLOCKII.ABSC[4] = Utilities::Text2Double(&ABSC4Box);
	agc.BLOCKII.SLOPE[0] = Utilities::Text2Double(&SLOPE0Box);
	agc.BLOCKII.SLOPE[1] = Utilities::Text2Double(&SLOPE1Box);
	agc.BLOCKII.SLOPE[2] = Utilities::Text2Double(&SLOPE2Box);
	agc.BLOCKII.SLOPE[3] = Utilities::Text2Double(&SLOPE3Box);
	agc.BLOCKII.SLOPE[4] = Utilities::Text2Double(&SLOPE4Box);
	agc.BLOCKII.IGNAOSQ = Utilities::Text2Double(&IGNAOSQBox);
	agc.BLOCKII.IGNAOSR = Utilities::Text2Double(&IGNAOSRBox);
	agc.BLOCKII.VIGN = Utilities::Text2Double(&VIGNBox);
	agc.BLOCKII.RIGNX = Utilities::Text2Double(&RIGNXBox);
	agc.BLOCKII.RIGNZ = Utilities::Text2Double(&RIGNZBox);
	agc.BLOCKII.KIGNXB4 = Utilities::Text2Double(&KIGNXBox);
	agc.BLOCKII.KIGNYB8 = Utilities::Text2Double(&KIGNYBox);
	agc.BLOCKII.KIGNVB4 = Utilities::Text2Double(&KIGNVBox);
	agc.BLOCKII.J1PARM = Utilities::Text2Double(&J1PARMBox);
	agc.BLOCKII.K1PARM = Utilities::Text2Double(&K1PARMBox);
	agc.BLOCKII.J2PARM = Utilities::Text2Double(&J2PARMBox);
	agc.BLOCKII.K2PARM = Utilities::Text2Double(&K2PARMBox);
	agc.BLOCKII.THETCRIT = Utilities::Text2Double(&THETCRITBox);
	agc.BLOCKII.RAMIN = Utilities::Text2Double(&RAMINBox);
	agc.BLOCKII.DELTTFAP = Utilities::Text2Double(&DELTTFAPBox);

	RopeNameBox.GetWindowText(string);
	std::wstring ws = std::wstring(string.GetString());
	agc.RopeName = std::string(ws.begin(), ws.end());

	PIOSDataSetBox.GetWindowText(string);
	ws = std::wstring(string.GetString());
	agc.PIOSDataSetName = std::string(ws.begin(), ws.end());

	agc.BLOCKII.R2Model = (R2ModelBox.GetCheck() != 0);

	agc.LGCDATA.IMUBiasCompensation.PBIASX = Utilities::Text2Double(&PBIASXBox);
	agc.LGCDATA.IMUBiasCompensation.PIPASCFX = Utilities::Text2Double(&PIPASCFXBox);
	agc.LGCDATA.IMUBiasCompensation.PBIASY = Utilities::Text2Double(&PBIASYBox);
	agc.LGCDATA.IMUBiasCompensation.PIPASCFY = Utilities::Text2Double(&PIPASCFYBox);
	agc.LGCDATA.IMUBiasCompensation.PBIASZ = Utilities::Text2Double(&PBIASZBox);
	agc.LGCDATA.IMUBiasCompensation.PIPASCFZ = Utilities::Text2Double(&PIPASCFZBox);
	agc.LGCDATA.IMUBiasCompensation.NBDX = Utilities::Text2Double(&NBDXBox);
	agc.LGCDATA.IMUBiasCompensation.NBDY = Utilities::Text2Double(&NBDYBox);
	agc.LGCDATA.IMUBiasCompensation.NBDZ = Utilities::Text2Double(&NBDZBox);
	agc.LGCDATA.IMUBiasCompensation.ADIAX = Utilities::Text2Double(&ADIAXBox);
	agc.LGCDATA.IMUBiasCompensation.ADIAY = Utilities::Text2Double(&ADIAYBox);
	agc.LGCDATA.IMUBiasCompensation.ADIAZ = Utilities::Text2Double(&ADIAZBox);
	agc.LGCDATA.IMUBiasCompensation.ADSRAX = Utilities::Text2Double(&ADSRAXBox);
	agc.LGCDATA.IMUBiasCompensation.ADSRAY = Utilities::Text2Double(&ADSRAYBox);
	agc.LGCDATA.IMUBiasCompensation.ADSRAZ = Utilities::Text2Double(&ADSRAZBox);

	int message = agc.RunLGC();

	//Write output
	switch (message)
	{
	case 0:
		OutputBox.SetWindowTextW(L"Padload generation successful!");
		break;
	case 1:
		OutputBox.SetWindowTextW(L"Rope not found!");
		break;
	case 2:
		OutputBox.SetWindowTextW(L"Rope not supported!");
		break;
	case 3:
		OutputBox.SetWindowTextW(L"LaunchMJD not supported by rope!");
		break;
	case 4:
		OutputBox.SetWindowTextW(L"PIOS data set not found!");
		break;
	}
}


void LGCPadloadGenerator::OnCbnSelchangeCombo1()
{
	UpdateRopeSpecificEditFields();
	UpdateData(FALSE);
}

void LGCPadloadGenerator::OnCbnSelchangeCombo2()
{
	MissionBox.GetLBText(MissionBox.GetCurSel(), MissionNameValue);
	UpdateData(FALSE);

	//Load some defaults
	DockedBox.SetCheck(BST_UNCHECKED);
	ROLLTIMEBox.SetWindowText(L"6.0");
	PITCHTIMEBox.SetWindowText(L"6.0");

	//Select new rope first
	switch (MissionBox.GetCurSel())
	{
	case 1: //Apollo 9
		RopeNameBox.SetCurSel(LGC_SUNDANCE306);
		break;
	case 2: //Apollo 10
		RopeNameBox.SetCurSel(LGC_LUMINARY069R2);
		break;
	case 3: //Apollo 11
		RopeNameBox.SetCurSel(LGC_LUMINARY099);
		break;
	case 4: //Apollo 12
		RopeNameBox.SetCurSel(LGC_LUMINARY116);
		break;
	case 5: //Apollo 13
		RopeNameBox.SetCurSel(LGC_LUMINARY131R1);
		break;
	case 6: //Apollo 14
		RopeNameBox.SetCurSel(LGC_LUMINARY178);
		break;
	case 7: //Apollo 15
	case 8: //Apollo 16
	case 9: //Apollo 17
		RopeNameBox.SetCurSel(LGC_LUMINARY210);
		break;
	}

	//Update to new rope here, so mission specific numbers can be overwritten later
	UpdateRopeSpecificEditFields();

	//Now mission specific numbers
	switch (MissionBox.GetCurSel())
	{
	case 1: //Apollo 9
		Apollo9Padload();
		break;
	case 2: //Apollo 10
		Apollo10Padload();
		break;
	case 3: //Apollo 11
		Apollo11Padload();
		break;
	case 4: //Apollo 12
		Apollo12Padload();
		break;
	case 5: //Apollo 13
		LaunchMJDInput.SetWindowText(L"40687.8006944444444");
		LSLatitudeBox.SetWindowText(L"-3.6686");
		LSLongitudeBox.SetWindowText(L"-17.4842");
		LSAltitudeBox.SetWindowText(L"-1405.0");
		LMMassBox.SetWindowTextW(L"33872.3");
		CSMMassBox.SetWindowTextW(L"37580.3");
		DockedBox.SetCheck(BST_UNCHECKED);
		HIASCENTBox.SetWindowText(L"10900.0");
		WRENDPOSBox.SetWindowTextW(L"10000.0");
		WRENDVELBox.SetWindowTextW(L"10.0");
		WSHAFTBox.SetWindowTextW(L"15.0");
		WTRUNBox.SetWindowTextW(L"15.0");
		RMAXBox.SetWindowTextW(L"2000.0");
		VMAXBox.SetWindowTextW(L"2.0");
		WSURFPOSBox.SetWindowText(L"0");
		WSURFVELBox.SetWindowText(L"0");
		SHAFTVARBox.SetWindowTextW(L"1.0");
		TRUNVARBox.SetWindowTextW(L"1.0");
		AGSKBox.SetWindowText(L"100.0");
		TLANDBox.SetWindowText(L"103.7433811");
		IGNAOSQBox.SetWindowText(L"7.63");
		IGNAOSRBox.SetWindowText(L"0.57");
		VIGNBox.SetWindowText(L"5545.3644");
		RIGNXBox.SetWindowText(L"-133371.54");
		RIGNZBox.SetWindowText(L"-1445069.5");
		KIGNXBox.SetWindowText(L"-0.331");
		KIGNYBox.SetWindowText(L"-5.8694e-7");
		KIGNVBox.SetWindowText(L"-438.0");
		J1PARMBox.SetWindowText(L"6042735.9");
		K1PARMBox.SetWindowText(L"-3.1743891e5");
		J2PARMBox.SetWindowText(L"6046910.4");
		K2PARMBox.SetWindowText(L"-6.2459985e5");
		THETCRITBox.SetWindowText(L"-17.183277");
		RAMINBox.SetWindowText(L"5.88048494e6");
		DELTTFAPBox.SetWindowText(L"-90.0");
		break;
	case 6: //Apollo 14
		LaunchMJDInput.SetWindowTextW(L"40982.84930555555");
		LSAltitudeBox.SetWindowTextW(L"-1405.2");
		LSLatitudeBox.SetWindowTextW(L"-3.67329493");
		LSLongitudeBox.SetWindowTextW(L"-17.46428902");
		LMMassBox.SetWindowTextW(L"34150.2");
		CSMMassBox.SetWindowTextW(L"36524.5");
		DockedBox.SetCheck(BST_UNCHECKED);
		HIASCENTBox.SetWindowText(L"10900.0");
		WRENDPOSBox.SetWindowTextW(L"10000.0");
		WRENDVELBox.SetWindowTextW(L"10.0");
		WSHAFTBox.SetWindowTextW(L"15.0");
		WTRUNBox.SetWindowTextW(L"15.0");
		RMAXBox.SetWindowTextW(L"2000.0");
		VMAXBox.SetWindowTextW(L"2.0");
		WSURFPOSBox.SetWindowText(L"0");
		WSURFVELBox.SetWindowText(L"0");
		SHAFTVARBox.SetWindowTextW(L"1.0");
		TRUNVARBox.SetWindowTextW(L"1.0");
		AGSKBox.SetWindowText(L"100.0");
		TLANDBox.SetWindowText(L"108.9202417");
		ABSC0Box.SetWindowText(L"-238000.0");
		ABSC1Box.SetWindowText(L"-57000.0");
		ABSC2Box.SetWindowText(L"-49000.0");
		ABSC3Box.SetWindowText(L"-11200.0");
		ABSC4Box.SetWindowText(L"-5000.0");
		SLOPE0Box.SetWindowText(L"-1.105e-2");
		SLOPE1Box.SetWindowText(L"-1.088e-1");
		SLOPE2Box.SetWindowText(L"3.704e-2");
		SLOPE3Box.SetWindowText(L"-7.903e-2");
		SLOPE4Box.SetWindowText(L"-1.2e-2");
		IGNAOSQBox.SetWindowText(L"6.027");
		IGNAOSRBox.SetWindowText(L"-0.016");
		VIGNBox.SetWindowText(L"5546.447179");
		RIGNXBox.SetWindowText(L"-140345.7283");
		RIGNZBox.SetWindowText(L"-1464979.987");
		KIGNXBox.SetWindowText(L"-0.419");
		KIGNYBox.SetWindowText(L"-9.05e-7");
		KIGNVBox.SetWindowText(L"-470.0");
		J1PARMBox.SetWindowText(L"6.0469527e6");
		K1PARMBox.SetWindowText(L"-3.1502779e5");
		J2PARMBox.SetWindowText(L"6.0486099e6");
		K2PARMBox.SetWindowText(L"-6.2763026e5");
		THETCRITBox.SetWindowText(L"-17.41421853");
		RAMINBox.SetWindowText(L"5.88006825e6");
		DELTTFAPBox.SetWindowText(L"-90.0");
		break;
	case 7: //Apollo 15
		LaunchMJDInput.SetWindowTextW(L"41158.565277778");
		LSAltitudeBox.SetWindowTextW(L"-3550.0");
		LSLatitudeBox.SetWindowTextW(L"26.073889");
		LSLongitudeBox.SetWindowTextW(L"3.6538889");
		LMMassBox.SetWindowTextW(L"36697.6");
		CSMMassBox.SetWindowTextW(L"37921.9");
		DockedBox.SetCheck(BST_UNCHECKED);
		HIASCENTBox.SetWindowText(L"10900.0");
		WRENDPOSBox.SetWindowTextW(L"10000.0");
		WRENDVELBox.SetWindowTextW(L"10.0");
		WSHAFTBox.SetWindowTextW(L"15.0");
		WTRUNBox.SetWindowTextW(L"15.0");
		RMAXBox.SetWindowTextW(L"2000.0");
		VMAXBox.SetWindowTextW(L"2.0");
		WSURFPOSBox.SetWindowText(L"0");
		WSURFVELBox.SetWindowText(L"0");
		SHAFTVARBox.SetWindowTextW(L"1.0");
		TRUNVARBox.SetWindowTextW(L"1.0");
		AGSKBox.SetWindowText(L"100.0");
		TLANDBox.SetWindowText(L"104.6824861");
		ABSC0Box.SetWindowText(L"-593000.0");
		ABSC1Box.SetWindowText(L"-196000.0");
		ABSC2Box.SetWindowText(L"-106000.0");
		ABSC3Box.SetWindowText(L"-89400.0");
		ABSC4Box.SetWindowText(L"-29200.0");
		SLOPE0Box.SetWindowText(L"2.0968e-2");
		SLOPE1Box.SetWindowText(L"0.0");
		SLOPE2Box.SetWindowText(L"1.86747e-1");
		SLOPE3Box.SetWindowText(L"-2.14784e-1");
		SLOPE4Box.SetWindowText(L"1.13014e-2");
		IGNAOSQBox.SetWindowText(L"5.961");
		IGNAOSRBox.SetWindowText(L"0.163");
		VIGNBox.SetWindowText(L"5548.14101");
		RIGNXBox.SetWindowText(L"-162539.6686");
		RIGNZBox.SetWindowText(L"-1547120.997");
		KIGNXBox.SetWindowText(L"-0.334");
		KIGNYBox.SetWindowText(L"-2.207e-7");
		KIGNVBox.SetWindowText(L"-498.0");
		J1PARMBox.SetWindowText(L"6.0410278e6");
		K1PARMBox.SetWindowText(L"-3.1137525e5");
		J2PARMBox.SetWindowText(L"6.0422821e6");
		K2PARMBox.SetWindowText(L"-6.2221362e5");
		THETCRITBox.SetWindowText(L"-17.510696");
		RAMINBox.SetWindowText(L"5.87303149e6");
		DELTTFAPBox.SetWindowText(L"-90.0");
		break;
	case 8: //Apollo 16
		LaunchMJDInput.SetWindowTextW(L"41423.7458333333");
		LSAltitudeBox.SetWindowTextW(L"-260.0");
		LSLatitudeBox.SetWindowTextW(L"-9.00028");
		LSLongitudeBox.SetWindowTextW(L"15.516389");
		LMMassBox.SetWindowTextW(L"36685.2");
		CSMMassBox.SetWindowTextW(L"39354.1");
		DockedBox.SetCheck(BST_UNCHECKED);
		HIASCENTBox.SetWindowText(L"10900.0");
		WRENDPOSBox.SetWindowTextW(L"10000.0");
		WRENDVELBox.SetWindowTextW(L"10.0");
		WSHAFTBox.SetWindowTextW(L"15.0");
		WTRUNBox.SetWindowTextW(L"15.0");
		RMAXBox.SetWindowTextW(L"2000.0");
		VMAXBox.SetWindowTextW(L"2.0");
		WSURFPOSBox.SetWindowText(L"0");
		WSURFVELBox.SetWindowText(L"0");
		SHAFTVARBox.SetWindowTextW(L"1.0");
		TRUNVARBox.SetWindowTextW(L"1.0");
		AGSKBox.SetWindowText(L"90.0");
		TLANDBox.SetWindowText(L"98.7784");
		ABSC0Box.SetWindowText(L"-692000.0");
		ABSC1Box.SetWindowText(L"-524000.0");
		ABSC2Box.SetWindowText(L"-234000.0");
		ABSC3Box.SetWindowText(L"-164000.0");
		ABSC4Box.SetWindowText(L"-100000.0");
		SLOPE0Box.SetWindowText(L"9.821428e-2");
		SLOPE1Box.SetWindowText(L"-5.17241e-3");
		SLOPE2Box.SetWindowText(L"-6.428571e-2");
		SLOPE3Box.SetWindowText(L"1.5625e-2");
		SLOPE4Box.SetWindowText(L"0.0");
		IGNAOSQBox.SetWindowText(L"5.442");
		IGNAOSRBox.SetWindowText(L"0.094");
		VIGNBox.SetWindowText(L"5543.4605");
		RIGNXBox.SetWindowText(L"-159548.72");
		RIGNZBox.SetWindowText(L"-1547623.3");
		KIGNXBox.SetWindowText(L"-0.334");
		KIGNYBox.SetWindowText(L"-2.207e-7");
		KIGNVBox.SetWindowText(L"-498.0");
		J1PARMBox.SetWindowText(L"6.051373e6");
		K1PARMBox.SetWindowText(L"-6.055552e5");
		J2PARMBox.SetWindowText(L"6.0383801e6");
		K2PARMBox.SetWindowText(L"-3.1696236e5");
		THETCRITBox.SetWindowText(L"5.63040073");
		RAMINBox.SetWindowText(L"5.883824672e6");
		DELTTFAPBox.SetWindowText(L"-70.0");
		break;
	case 9: //Apollo 17
		LaunchMJDInput.SetWindowTextW(L"41658.120138888");
		LSAltitudeBox.SetWindowTextW(L"-3606.0");
		LSLatitudeBox.SetWindowTextW(L"20.164029");
		LSLongitudeBox.SetWindowTextW(L"30.749532");
		LMMassBox.SetWindowTextW(L"36759.3");
		CSMMassBox.SetWindowTextW(L"38115.5");
		DockedBox.SetCheck(BST_UNCHECKED);
		HIASCENTBox.SetWindowText(L"10900.0");
		WRENDPOSBox.SetWindowTextW(L"10000.0");
		WRENDVELBox.SetWindowTextW(L"10.0");
		WSHAFTBox.SetWindowTextW(L"15.0");
		WTRUNBox.SetWindowTextW(L"15.0");
		RMAXBox.SetWindowTextW(L"2000.0");
		VMAXBox.SetWindowTextW(L"2.0");
		WSURFPOSBox.SetWindowText(L"0");
		WSURFVELBox.SetWindowText(L"0");
		SHAFTVARBox.SetWindowTextW(L"1.0");
		TRUNVARBox.SetWindowTextW(L"1.0");
		AGSKBox.SetWindowText(L"110.0");
		TLANDBox.SetWindowText(L"113.0272472");
		ABSC0Box.SetWindowText(L"-205500.0");
		ABSC1Box.SetWindowText(L"-57500.0");
		ABSC2Box.SetWindowText(L"-44000.0");
		ABSC3Box.SetWindowText(L"-29000.0");
		ABSC4Box.SetWindowText(L"-16300.0");
		SLOPE0Box.SetWindowText(L"3.37837e-3");
		SLOPE1Box.SetWindowText(L"-3.40741e-1");
		SLOPE2Box.SetWindowText(L"1.333333e-1");
		SLOPE3Box.SetWindowText(L"-2.07086e-1");
		SLOPE4Box.SetWindowText(L"1.411e-2");
		IGNAOSQBox.SetWindowText(L"5.335");
		IGNAOSRBox.SetWindowText(L"0.003");
		VIGNBox.SetWindowText(L"5542.8976");
		RIGNXBox.SetWindowText(L"-156145.03");
		RIGNZBox.SetWindowText(L"-1541941.8");
		KIGNXBox.SetWindowText(L"-0.334");
		KIGNYBox.SetWindowText(L"-2.207e-7");
		KIGNVBox.SetWindowText(L"-498.0");
		J1PARMBox.SetWindowText(L"6.0457376e6");
		K1PARMBox.SetWindowText(L"-6.0598187e5");
		J2PARMBox.SetWindowText(L"6.036605e6");
		K2PARMBox.SetWindowText(L"-3.180295e5");
		THETCRITBox.SetWindowText(L"6.432346783");
		RAMINBox.SetWindowText(L"5.872844816e6");
		DELTTFAPBox.SetWindowText(L"-70.0");

		/*
		PBIASXBox.SetWindowText(L"1.64");		// cm/sec^2
		PIPASCFXBox.SetWindowText(L"-980.0");	// ppm
		PBIASYBox.SetWindowText(L"1.73");		// cm/sec^2
		PIPASCFYBox.SetWindowText(L"-560.0");	// ppm
		PBIASZBox.SetWindowText(L"1.6");		// cm/sec^2
		PIPASCFZBox.SetWindowText(L"-460.0");	// ppm
		NBDXBox.SetWindowText(L"0.1");			// meru
		NBDYBox.SetWindowText(L"0.4");			// meru
		NBDZBox.SetWindowText(L"-1.1");			// meru
		ADIAXBox.SetWindowText(L"12.0");		// meru/g
		ADIAYBox.SetWindowText(L"-4.0");		// meru/g
		ADIAZBox.SetWindowText(L"4.0");			// meru/g
		ADSRAXBox.SetWindowText(L"4.0");		// meru/g
		ADSRAYBox.SetWindowText(L"6.0");		// meru/g
		ADSRAZBox.SetWindowText(L"-8.0");		// meru/g
		*/
		break;
	}

	UpdateTotalMass();
}


void LGCPadloadGenerator::OnBnClickedCancel()
{
	CDialogEx::OnCancel();
}


void LGCPadloadGenerator::OnEnChangeEdit6()
{
	UpdateTotalMass();
}

void LGCPadloadGenerator::OnEnChangeEdit7()
{
	UpdateTotalMass();
}

void LGCPadloadGenerator::OnBnClickedCheck1()
{
	UpdateTotalMass();
}

void LGCPadloadGenerator::UpdateTotalMass()
{
	double TotalMass = 0.0;

	TotalMass += Utilities::Text2Double(&LMMassBox);

	if (DockedBox.GetCheck())
	{
		TotalMass += Utilities::Text2Double(&CSMMassBox);
	}

	Utilities::Double2Text(TotalMass, &TotalMassBox, 5);
}

void LGCPadloadGenerator::UpdateRopeSpecificEditFields()
{
	if (RopeNameBox.GetCurSel() == LGC_SUNDANCE306)
	{
		WSURFPOSBox.SetReadOnly(true);
		WSURFVELBox.SetReadOnly(true);
		VIGNBox.SetReadOnly(true);
		RIGNXBox.SetReadOnly(true);
		RIGNZBox.SetReadOnly(true);
		KIGNXBox.SetReadOnly(true);
		KIGNYBox.SetReadOnly(true);
		KIGNVBox.SetReadOnly(true);
	}
	else
	{
		WSURFPOSBox.SetReadOnly(false);
		WSURFVELBox.SetReadOnly(false);
		VIGNBox.SetReadOnly(false);
		RIGNXBox.SetReadOnly(false);
		RIGNZBox.SetReadOnly(false);
		KIGNXBox.SetReadOnly(false);
		KIGNYBox.SetReadOnly(false);
		KIGNVBox.SetReadOnly(false);
	}

	if (RopeNameBox.GetCurSel() == LGC_LUMINARY178 || RopeNameBox.GetCurSel() == LGC_LUMINARY210 || RopeNameBox.GetCurSel() == LGC_ZERLINA056)
	{
		ABSC0Box.SetReadOnly(false); ABSC1Box.SetReadOnly(false); ABSC2Box.SetReadOnly(false); ABSC3Box.SetReadOnly(false); ABSC4Box.SetReadOnly(false);
		SLOPE0Box.SetReadOnly(false); SLOPE1Box.SetReadOnly(false); SLOPE2Box.SetReadOnly(false); SLOPE3Box.SetReadOnly(false); SLOPE4Box.SetReadOnly(false);
	}
	else
	{
		ABSC0Box.SetReadOnly(true); ABSC1Box.SetReadOnly(true); ABSC2Box.SetReadOnly(true); ABSC3Box.SetReadOnly(true); ABSC4Box.SetReadOnly(true);
		SLOPE0Box.SetReadOnly(true); SLOPE1Box.SetReadOnly(true); SLOPE2Box.SetReadOnly(true); SLOPE3Box.SetReadOnly(true); SLOPE4Box.SetReadOnly(true);
	}

	if (RopeNameBox.GetCurSel() <= LGC_LUMINARY099)
	{
		J1PARMBox.SetReadOnly(true);
		K1PARMBox.SetReadOnly(true);
		J2PARMBox.SetReadOnly(true);
		K2PARMBox.SetReadOnly(true);
		THETCRITBox.SetReadOnly(true);
		RAMINBox.SetReadOnly(true);
	}
	else
	{
		J1PARMBox.SetReadOnly(false);
		K1PARMBox.SetReadOnly(false);
		J2PARMBox.SetReadOnly(false);
		K2PARMBox.SetReadOnly(false);
		THETCRITBox.SetReadOnly(false);
		RAMINBox.SetReadOnly(false);
	}

	//PIOS Data Set
	switch (RopeNameBox.GetCurSel())
	{
	case LGC_SUNDANCE306:
	case LGC_LUMINARY069:
	case LGC_LUMINARY069R2:
		PIOSDataSetBox.SetWindowTextW(L"NBY1969");
		break;
	case LGC_LUMINARY099:
		PIOSDataSetBox.SetWindowTextW(L"NBY1970_V1");
		break;
	case LGC_LUMINARY116:
	case LGC_LUMINARY131R1:
		PIOSDataSetBox.SetWindowTextW(L"NBY1970_V2");
		break;
	case LGC_LUMINARY178:
	case LGC_ZERLINA056:
		PIOSDataSetBox.SetWindowTextW(L"NBY1971");
		break;
	default:
		PIOSDataSetBox.SetWindowTextW(L"NBY1972");
		break;
	}
}

void LGCPadloadGenerator::Apollo9Padload()
{
	LaunchMJDInput.SetWindowTextW(L"40283.6666667");
	LMMassBox.SetWindowTextW(L"32401.2");
	CSMMassBox.SetWindowTextW(L"30052.6");
	HIASCENTBox.SetWindowTextW(L"10145.4");
	DockedBox.SetCheck(BST_CHECKED);
	WRENDPOSBox.SetWindowTextW(L"1000.0");
	WRENDVELBox.SetWindowTextW(L"1.0");
	WSHAFTBox.SetWindowTextW(L"5.0");
	WTRUNBox.SetWindowTextW(L"5.0");
	RMAXBox.SetWindowTextW(L"5000.0");
	VMAXBox.SetWindowTextW(L"5.0");
	SHAFTVARBox.SetWindowTextW(L"1.0");
	TRUNVARBox.SetWindowTextW(L"1.0");
	AGSKBox.SetWindowTextW(L"40.0");
	ROLLTIMEBox.SetWindowTextW(L"6.95");
	PITCHTIMEBox.SetWindowTextW(L"6.43");

	PBIASXBox.SetWindowText(L"0.31");		// cm/sec^2
	PIPASCFXBox.SetWindowText(L"-968.0");	// ppm
	PBIASYBox.SetWindowText(L"0.19");		// cm/sec^2
	PIPASCFYBox.SetWindowText(L"-941.0");	// ppm
	PBIASZBox.SetWindowText(L"0.0");		// cm/sec^2
	PIPASCFZBox.SetWindowText(L"-852.0");	// ppm
	NBDXBox.SetWindowText(L"4.6");			// meru
	NBDYBox.SetWindowText(L"5.0");			// meru
	NBDZBox.SetWindowText(L"4.5");			// meru
	ADIAXBox.SetWindowText(L"5.4");			// meru/g
	ADIAYBox.SetWindowText(L"-0.3");		// meru/g
	ADIAZBox.SetWindowText(L"19.6");		// meru/g
	ADSRAXBox.SetWindowText(L"-0.5");		// meru/g
	ADSRAYBox.SetWindowText(L"16.3");		// meru/g
	ADSRAZBox.SetWindowText(L"-1.7");		// meru/g
}

void LGCPadloadGenerator::Apollo10Padload()
{
	LaunchMJDInput.SetWindowTextW(L"40359.7006944");
	LSAltitudeBox.SetWindowTextW(L"-3073.26");
	LSLatitudeBox.SetWindowTextW(L"0.71388");
	LSLongitudeBox.SetWindowTextW(L"23.707773");
	LMMassBox.SetWindowTextW(L"31132.3");
	CSMMassBox.SetWindowTextW(L"37605.0");
	TotalMassBox.SetWindowTextW(L"31132.3");
	DockedBox.SetCheck(BST_UNCHECKED);
	HIASCENTBox.SetWindowTextW(L"8422.9");
	WRENDPOSBox.SetWindowTextW(L"10000.0");
	WRENDVELBox.SetWindowTextW(L"10.0");
	WSHAFTBox.SetWindowTextW(L"15.0");
	WTRUNBox.SetWindowTextW(L"15.0");
	RMAXBox.SetWindowTextW(L"2000.0");
	VMAXBox.SetWindowTextW(L"2.0");
	WSURFPOSBox.SetWindowText(L"52493.4383202");
	WSURFVELBox.SetWindowText(L"9.842519685");
	SHAFTVARBox.SetWindowTextW(L"1.0");
	TRUNVARBox.SetWindowTextW(L"1.0");
	AGSKBox.SetWindowTextW(L"90.0");
	TLANDBox.SetWindowText(L"100.847");
	VIGNBox.SetWindowText(L"5545.46");
	RIGNXBox.SetWindowText(L"-130519.86");
	RIGNZBox.SetWindowText(L"-1430097.4");
	KIGNXBox.SetWindowText(L"-0.617631");
	KIGNYBox.SetWindowText(L"-0.755e-6");
	KIGNVBox.SetWindowText(L"-410.0");

	PBIASXBox.SetWindowText(L"-0.41");		// cm/sec^2
	PIPASCFXBox.SetWindowText(L"-430.0");	// ppm
	PBIASYBox.SetWindowText(L"0.18");		// cm/sec^2
	PIPASCFYBox.SetWindowText(L"-840.0");	// ppm
	PBIASZBox.SetWindowText(L"-0.03");		// cm/sec^2
	PIPASCFZBox.SetWindowText(L"-530.0");	// ppm
	NBDXBox.SetWindowText(L"-3.2");			// meru
	NBDYBox.SetWindowText(L"1.5");			// meru
	NBDZBox.SetWindowText(L"-1.2");			// meru
	ADIAXBox.SetWindowText(L"1.0");			// meru/g
	ADIAYBox.SetWindowText(L"20.0");		// meru/g
	ADIAZBox.SetWindowText(L"-24.0");		// meru/g
	ADSRAXBox.SetWindowText(L"5.0");		// meru/g
	ADSRAYBox.SetWindowText(L"2.0");		// meru/g
	ADSRAZBox.SetWindowText(L"-1.0");		// meru/g
}

void LGCPadloadGenerator::Apollo11Padload()
{
	LaunchMJDInput.SetWindowTextW(L"40418.5638889");
	LSAltitudeBox.SetWindowTextW(L"-2671.9684");
	LSLatitudeBox.SetWindowTextW(L"0.691395");
	LSLongitudeBox.SetWindowTextW(L"23.71689");
	LMMassBox.SetWindowTextW(L"33663.9535");
	CSMMassBox.SetWindowTextW(L"36470.4");
	TotalMassBox.SetWindowTextW(L"33663.9535");
	DockedBox.SetCheck(BST_UNCHECKED);
	HIASCENTBox.SetWindowTextW(L"11075.1");
	WRENDPOSBox.SetWindowTextW(L"10000.0");
	WRENDVELBox.SetWindowTextW(L"10.0");
	WSHAFTBox.SetWindowTextW(L"15.0");
	WTRUNBox.SetWindowTextW(L"15.0");
	RMAXBox.SetWindowTextW(L"2000.0");
	VMAXBox.SetWindowTextW(L"2.0");
	WSURFPOSBox.SetWindowText(L"5000");
	WSURFVELBox.SetWindowText(L"5");
	SHAFTVARBox.SetWindowTextW(L"1.0");
	TRUNVARBox.SetWindowTextW(L"1.0");
	AGSKBox.SetWindowTextW(L"90.0");
	TLANDBox.SetWindowText(L"102.786402777777778");
	IGNAOSQBox.SetWindowText(L"6.25");
	IGNAOSRBox.SetWindowText(L"0.63");
	VIGNBox.SetWindowText(L"5545.46");
	RIGNXBox.SetWindowText(L"-130519.86");
	RIGNZBox.SetWindowText(L"-1432597.3");
	KIGNXBox.SetWindowText(L"-0.617631");
	KIGNYBox.SetWindowText(L"-0.755e-6");
	KIGNVBox.SetWindowText(L"-410.0");
	DELTTFAPBox.SetWindowText(L"-110.0");

	PBIASXBox.SetWindowText(L"-0.12");		// cm/sec^2
	PIPASCFXBox.SetWindowText(L"-210.0");	// ppm
	PBIASYBox.SetWindowText(L"0.17");		// cm/sec^2
	PIPASCFYBox.SetWindowText(L"-210.0");	// ppm
	PBIASZBox.SetWindowText(L"0.07");		// cm/sec^2
	PIPASCFZBox.SetWindowText(L"-920.0");	// ppm
	NBDXBox.SetWindowText(L"1.3");			// meru
	NBDYBox.SetWindowText(L"-1.5");			// meru
	NBDZBox.SetWindowText(L"-1.9");			// meru
	ADIAXBox.SetWindowText(L"57.0");		// meru/g
	ADIAYBox.SetWindowText(L"-4.0");		// meru/g
	ADIAZBox.SetWindowText(L"20.0");		// meru/g
	ADSRAXBox.SetWindowText(L"2.0");		// meru/g
	ADSRAYBox.SetWindowText(L"26.0");		// meru/g
	ADSRAZBox.SetWindowText(L"-11.0");		// meru/g
}

void LGCPadloadGenerator::Apollo12Padload()
{
	LaunchMJDInput.SetWindowTextW(L"40539.6819444");
	LSAltitudeBox.SetWindowTextW(L"-2371.27");
	LSLatitudeBox.SetWindowTextW(L"-2.9822165");
	LSLongitudeBox.SetWindowTextW(L"-23.391933");
	LMMassBox.SetWindowTextW(L"33928.8");
	CSMMassBox.SetWindowTextW(L"37232.2");
	DockedBox.SetCheck(BST_UNCHECKED);
	HIASCENTBox.SetWindowText(L"10900.0");
	WRENDPOSBox.SetWindowTextW(L"10000.0");
	WRENDVELBox.SetWindowTextW(L"10.0");
	WSHAFTBox.SetWindowTextW(L"15.0");
	WTRUNBox.SetWindowTextW(L"15.0");
	RMAXBox.SetWindowTextW(L"2000.0");
	VMAXBox.SetWindowTextW(L"2.0");
	WSURFPOSBox.SetWindowText(L"0");
	WSURFVELBox.SetWindowText(L"0");
	SHAFTVARBox.SetWindowTextW(L"1.0");
	TRUNVARBox.SetWindowTextW(L"1.0");
	AGSKBox.SetWindowTextW(L"100.0");
	TLANDBox.SetWindowText(L"110.57895");
	IGNAOSQBox.SetWindowText(L"7.97");
	IGNAOSRBox.SetWindowText(L"0.56");
	VIGNBox.SetWindowText(L"5551.1299");
	RIGNXBox.SetWindowText(L"-133067.52");
	RIGNZBox.SetWindowText(L"-1437887.4");
	KIGNXBox.SetWindowText(L"-0.331");
	KIGNYBox.SetWindowText(L"-5.8694e-7");
	KIGNVBox.SetWindowText(L"-438.0");
	J1PARMBox.SetWindowText(L"6032567.5");
	K1PARMBox.SetWindowText(L"-6.2726125e5");
	J2PARMBox.SetWindowText(L"6030470.0");
	K2PARMBox.SetWindowText(L"-3.1835146e5");
	THETCRITBox.SetWindowText(L"8.384852304");
	RAMINBox.SetWindowText(L"5.8768997e6");
	DELTTFAPBox.SetWindowText(L"-90.0");

	PBIASXBox.SetWindowText(L"-0.38");		// cm/sec^2
	PIPASCFXBox.SetWindowText(L"-660.0");	// ppm
	PBIASYBox.SetWindowText(L"0.02");		// cm/sec^2
	PIPASCFYBox.SetWindowText(L"-720.0");	// ppm
	PBIASZBox.SetWindowText(L"0.62");		// cm/sec^2
	PIPASCFZBox.SetWindowText(L"-890.0");	// ppm
	NBDXBox.SetWindowText(L"0.1");			// meru
	NBDYBox.SetWindowText(L"0.8");			// meru
	NBDZBox.SetWindowText(L"3.0");			// meru
	ADIAXBox.SetWindowText(L"17.0");		// meru/g
	ADIAYBox.SetWindowText(L"-15.0");		// meru/g
	ADIAZBox.SetWindowText(L"13.0");		// meru/g
	ADSRAXBox.SetWindowText(L"-2.0");		// meru/g
	ADSRAYBox.SetWindowText(L"4.0");		// meru/g
	ADSRAZBox.SetWindowText(L"-2.0");		// meru/g
}