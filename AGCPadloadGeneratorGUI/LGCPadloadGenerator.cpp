// LGCPadloadGenerator.cpp
//

#include "pch.h"
#include "AGCPadloadGeneratorGUI.h"
#include "LGCPadloadGenerator.h"
#include "afxdialogex.h"


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
}


BEGIN_MESSAGE_MAP(LGCPadloadGenerator, CDialogEx)
	ON_BN_CLICKED(IDOK, &LGCPadloadGenerator::OnBnClickedOk)
	ON_CBN_SELCHANGE(IDC_COMBO1, &LGCPadloadGenerator::OnCbnSelchangeCombo1)
	ON_CBN_SELCHANGE(IDC_COMBO2, &LGCPadloadGenerator::OnCbnSelchangeCombo2)
	ON_BN_CLICKED(IDCANCEL, &LGCPadloadGenerator::OnBnClickedCancel)
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

	RopeNameBox.AddString(L"Luminary131R1");
	RopeNameBox.SetCurSel(0);

	T504LMBox.SetWindowText(L"4.5");
	UNITWBox.SetWindowText(L"-0.5");
	LSLatitudeBox.SetWindowText(L"0.0");
	LSLongitudeBox.SetWindowText(L"0.0");
	LSAltitudeBox.SetWindowText(L"0.0");
	EphemerisSpanBox.SetWindowText(L"14.5");
	LaunchMJDInput.SetWindowText(L"0.0");

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

	RopeNameBox.GetWindowText(string);
	std::wstring ws = std::wstring(string.GetString());
	agc.RopeName = std::string(ws.begin(), ws.end());

	agc.RunLGC();
}


void LGCPadloadGenerator::OnCbnSelchangeCombo1()
{
	RopeNameBox.GetLBText(RopeNameBox.GetCurSel(), RopeNameValue);
	UpdateData(FALSE);
}


void LGCPadloadGenerator::OnCbnSelchangeCombo2()
{
	MissionBox.GetLBText(MissionBox.GetCurSel(), MissionNameValue);
	UpdateData(FALSE);

	switch (MissionBox.GetCurSel())
	{
	case 1: //Apollo 9
		break;
	case 2: //Apollo 10
		break;
	case 3: //Apollo 11
		LaunchMJDInput.SetWindowTextW(L"40418.5638889");
		LSAltitudeBox.SetWindowTextW(L"-2671.9684");
		LSLatitudeBox.SetWindowTextW(L"0.691395");
		LSLongitudeBox.SetWindowTextW(L"23.71689");
		break;
	case 4: //Apollo 12
		LaunchMJDInput.SetWindowTextW(L"40539.6819444");
		LSAltitudeBox.SetWindowTextW(L"-2371.27");
		LSLatitudeBox.SetWindowTextW(L"-2.9822165");
		LSLongitudeBox.SetWindowTextW(L"-23.391933");
		break;
	case 5: //Apollo 13
		LaunchMJDInput.SetWindowText(L"40687.8006944444444");
		LSLatitudeBox.SetWindowText(L"-3.6686");
		LSLongitudeBox.SetWindowText(L"-17.4842");
		LSAltitudeBox.SetWindowText(L"-1405.0");
		break;
	}
}


void LGCPadloadGenerator::OnBnClickedCancel()
{
	CDialogEx::OnCancel();
}
