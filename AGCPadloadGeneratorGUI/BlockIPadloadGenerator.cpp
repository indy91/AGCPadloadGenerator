// BlockIPadloadGenerator.cpp : implementation file
//

#include "pch.h"
#include "AGCPadloadGeneratorGUI.h"
#include "BlockIPadloadGenerator.h"
#include "afxdialogex.h"


// BlockIPadloadGenerator dialog

IMPLEMENT_DYNAMIC(BlockIPadloadGenerator, CDialogEx)

BlockIPadloadGenerator::BlockIPadloadGenerator(CWnd* pParent /*=nullptr*/)
	: CDialogEx(IDD_DIALOG3, pParent)
{

}

BlockIPadloadGenerator::~BlockIPadloadGenerator()
{
}

void BlockIPadloadGenerator::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_COMBO1, Launchpad);
	DDX_Control(pDX, IDC_COMBO3, MissionBox);
	DDX_Control(pDX, IDC_COMBO2, RopeNameBox);
	DDX_Control(pDX, IDC_EDIT1, LaunchMJDInput);
	DDX_Control(pDX, IDC_EDIT2, LaunchAzimuthBox);
	DDX_Control(pDX, IDC_EDIT4, SPS1EccentricityBox);
	DDX_Control(pDX, IDC_EDIT7, SPS2EccentricityBox);
	DDX_Control(pDX, IDC_EDIT5, SPS1SMABox);
	DDX_Control(pDX, IDC_EDIT8, SPS2SMABox);
}


BEGIN_MESSAGE_MAP(BlockIPadloadGenerator, CDialogEx)
	ON_BN_CLICKED(IDOK, &BlockIPadloadGenerator::OnBnClickedOk)
	ON_CBN_SELCHANGE(IDC_COMBO3, &BlockIPadloadGenerator::OnCbnSelchangeCombo3)
END_MESSAGE_MAP()

BOOL BlockIPadloadGenerator::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	Launchpad.AddString(L"LC-34");
	Launchpad.AddString(L"LC-39A");
	Launchpad.AddString(L"LC-39B");
	Launchpad.SetCurSel(1);

	MissionBox.AddString(L"Manual");
	MissionBox.AddString(L"AS-202");
	MissionBox.AddString(L"Apollo 4");
	MissionBox.AddString(L"Apollo 6");
	MissionBox.SetCurSel(0);

	RopeNameBox.AddString(L"Corona261");
	RopeNameBox.AddString(L"Solarium055");
	RopeNameBox.SetCurSel(1);

	LaunchAzimuthBox.SetWindowText(L"72.0");

	return TRUE;
}

// BlockIPadloadGenerator message handlers


void BlockIPadloadGenerator::OnBnClickedOk()
{
	CString string;

	agc.LaunchMJD = Utilities::Text2Double(&LaunchMJDInput);
	agc.LaunchAzimuth = Utilities::Text2Double(&LaunchAzimuthBox);
	agc.BLOCKI.e_SPS1 = Utilities::Text2Double(&SPS1EccentricityBox);
	agc.BLOCKI.a_SPS1 = Utilities::Text2Double(&SPS1SMABox);
	agc.BLOCKI.e_SPS2 = Utilities::Text2Double(&SPS2EccentricityBox);
	agc.BLOCKI.a_SPS2 = Utilities::Text2Double(&SPS2SMABox);

	Launchpad.GetWindowText(string);
	std::wstring ws = std::wstring(string.GetString());
	agc.Pad = std::string(ws.begin(), ws.end());

	RopeNameBox.GetWindowText(string);
	ws = std::wstring(string.GetString());
	agc.RopeName = std::string(ws.begin(), ws.end());

	agc.RunBlockI();
}


void BlockIPadloadGenerator::OnCbnSelchangeCombo3()
{
	switch (MissionBox.GetCurSel())
	{
	case 1: //AS-202
		RopeNameBox.SetCurSel(0); //Corona 261
		LaunchMJDInput.SetWindowTextW(L"39362.719120");
		LaunchAzimuthBox.SetWindowTextW(L"105.0");
		SPS1EccentricityBox.SetWindowTextW(L"0.10988556");
		SPS2EccentricityBox.SetWindowTextW(L"0.25341222");
		SPS1SMABox.SetWindowTextW(L"6855280.0");
		SPS2SMABox.SetWindowTextW(L"8623082.5");
		break;
	case 2: //Apollo 4
		RopeNameBox.SetCurSel(1); //Solarium 55
		LaunchMJDInput.SetWindowTextW(L"39803.5");
		LaunchAzimuthBox.SetWindowTextW(L"72.0");
		SPS1EccentricityBox.SetWindowTextW(L"0.5934490037");
		SPS2EccentricityBox.SetWindowTextW(L"0.999071629063702");
		SPS1SMABox.SetWindowTextW(L"15487553.0");
		SPS2SMABox.SetWindowTextW(L"6891085630.5");
		break;
	case 3: //Apollo 6
		RopeNameBox.SetCurSel(1); //Solarium 55
		LaunchMJDInput.SetWindowTextW(L"39950.5");
		LaunchAzimuthBox.SetWindowTextW(L"72.0");
		SPS1EccentricityBox.SetWindowTextW(L"0.63429326");
		SPS2EccentricityBox.SetWindowTextW(L"1.0164386");
		SPS1SMABox.SetWindowTextW(L"17512783.0");
		SPS2SMABox.SetWindowTextW(L"-390049424.0");
		break;
	}
}
