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

BOOL LGCPadloadGenerator::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	MissionBox.AddString(L"Manual");
	MissionBox.AddString(L"Apollo 9");
	MissionBox.AddString(L"Apollo 10");
	MissionBox.AddString(L"Apollo 11");
	MissionBox.AddString(L"Apollo 12");
	MissionBox.AddString(L"Apollo 13");
	MissionBox.SetCurSel(0);

	RopeNameBox.AddString(L"Sundance306");
	RopeNameBox.AddString(L"Luminary069");
	RopeNameBox.AddString(L"Luminary069R2");
	RopeNameBox.AddString(L"Luminary099");
	RopeNameBox.AddString(L"Luminary116");
	RopeNameBox.AddString(L"Luminary131R1");
	RopeNameBox.SetCurSel(LGC_LUMINARY116);

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

	RopeNameBox.GetWindowText(string);
	std::wstring ws = std::wstring(string.GetString());
	agc.RopeName = std::string(ws.begin(), ws.end());

	agc.RunLGC();
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

	DockedBox.SetCheck(BST_UNCHECKED);
	ROLLTIMEBox.SetWindowText(L"6.0");
	PITCHTIMEBox.SetWindowText(L"6.0");

	switch (MissionBox.GetCurSel())
	{
	case 1: //Apollo 9
		RopeNameBox.SetCurSel(LGC_SUNDANCE306);
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
		break;
	case 2: //Apollo 10
		RopeNameBox.SetCurSel(LGC_LUMINARY069R2);
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
		break;
	case 3: //Apollo 11
		RopeNameBox.SetCurSel(LGC_LUMINARY099);
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
		break;
	case 4: //Apollo 12
		RopeNameBox.SetCurSel(LGC_LUMINARY116);
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
		break;
	case 5: //Apollo 13
		RopeNameBox.SetCurSel(LGC_LUMINARY131R1);
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
		break;
	}

	UpdateTotalMass();
	UpdateRopeSpecificEditFields();
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
	}
	else
	{
		WSURFPOSBox.SetReadOnly(false);
		WSURFVELBox.SetReadOnly(false);
	}
}