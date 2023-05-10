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
	DDX_Control(pDX, IDC_EDIT3, TAtlanticBox);
	DDX_Control(pDX, IDC_EDIT10, TPacificBox);
	DDX_Control(pDX, IDC_EDIT6, AtlanticLatitudeBox);
	DDX_Control(pDX, IDC_EDIT9, AtlanticLongitudeBox);
	DDX_Control(pDX, IDC_EDIT11, PacificLatitudeBox);
	DDX_Control(pDX, IDC_EDIT12, PacificLongitudeBox);
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
	agc.BLOCKI.T_ATL = Utilities::Text2Double(&TAtlanticBox);
	agc.BLOCKI.lat_ATL = Utilities::Text2Double(&AtlanticLatitudeBox);
	agc.BLOCKI.lng_ATL = Utilities::Text2Double(&AtlanticLongitudeBox);
	agc.BLOCKI.T_PAC = Utilities::Text2Double(&TPacificBox);
	agc.BLOCKI.lat_PAC = Utilities::Text2Double(&PacificLatitudeBox);
	agc.BLOCKI.lng_PAC = Utilities::Text2Double(&PacificLongitudeBox);

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
		Launchpad.SetCurSel(0); //LC-34
		LaunchMJDInput.SetWindowTextW(L"39362.719120");
		LaunchAzimuthBox.SetWindowTextW(L"105.0");
		SPS1EccentricityBox.SetWindowTextW(L"0.10988556");
		SPS2EccentricityBox.SetWindowTextW(L"0.25341222");
		SPS1SMABox.SetWindowTextW(L"6855280.0");
		SPS2SMABox.SetWindowTextW(L"8623082.5");
		TAtlanticBox.SetWindowTextW(L"1400.0");
		AtlanticLatitudeBox.SetWindowTextW(L"4.0");
		AtlanticLongitudeBox.SetWindowTextW(L"-31.0");
		TPacificBox.SetWindowTextW(L"5648.0");
		PacificLatitudeBox.SetWindowTextW(L"17.25");
		PacificLongitudeBox.SetWindowTextW(L"170.0");

		//Temporary, until there is input boxes
		agc.BLOCKI.TROLL = 8.0;
		agc.BLOCKI.TPITCH = 10.0;
		agc.BLOCKI.TENDPITCH = 126.0;

		//Measured from local vertical
		agc.BLOCKI.POLYCOFF[0] = 8.99456725e1 - 90.0 - 33.0; //Measured from zero, offset by 33° for the IMU
		agc.BLOCKI.POLYCOFF[1] = 6.05731859e-4;
		agc.BLOCKI.POLYCOFF[2] = -3.33462947e-3;
		agc.BLOCKI.POLYCOFF[3] = -1.8166406e-4;
		agc.BLOCKI.POLYCOFF[4] = 3.17822761e-6;
		agc.BLOCKI.POLYCOFF[5] = -1.88355082e-8;
		agc.BLOCKI.POLYCOFF[6] = 3.93873259e-11;

		break;
	case 2: //Apollo 4
		RopeNameBox.SetCurSel(1); //Solarium 55
		Launchpad.SetCurSel(1); //LC-39A
		LaunchMJDInput.SetWindowTextW(L"39803.5");
		LaunchAzimuthBox.SetWindowTextW(L"72.0");
		SPS1EccentricityBox.SetWindowTextW(L"0.5934490037");
		SPS2EccentricityBox.SetWindowTextW(L"0.999071629063702");
		SPS1SMABox.SetWindowTextW(L"15487553.0");
		SPS2SMABox.SetWindowTextW(L"6891085630.5");
		TAtlanticBox.SetWindowTextW(L"1400.0");
		TPacificBox.SetWindowTextW(L"30921.42");
		AtlanticLatitudeBox.SetWindowTextW(L"28.29028886");
		AtlanticLongitudeBox.SetWindowTextW(L"-19.5");
		PacificLatitudeBox.SetWindowTextW(L"30.04649677");
		PacificLongitudeBox.SetWindowTextW(L"-171.0");

		//Temporary, until there is input boxes
		agc.BLOCKI.POLYCOFF[0] = 0.4788289;
		agc.BLOCKI.POLYCOFF[1] = 0.3411991e-1;
		agc.BLOCKI.POLYCOFF[2] = 0.1178593e-1;
		agc.BLOCKI.POLYCOFF[3] = 0.5318752e-4;
		agc.BLOCKI.POLYCOFF[4] = -0.2937101e-5;
		agc.BLOCKI.POLYCOFF[5] = 0.232939e-7;
		agc.BLOCKI.POLYCOFF[6] = -0.5793149e-10;
		break;
	case 3: //Apollo 6
		RopeNameBox.SetCurSel(1); //Solarium 55
		Launchpad.SetCurSel(1); //LC-39A
		LaunchMJDInput.SetWindowTextW(L"39950.5");
		LaunchAzimuthBox.SetWindowTextW(L"72.0");
		SPS1EccentricityBox.SetWindowTextW(L"0.63429326");
		SPS2EccentricityBox.SetWindowTextW(L"1.0164386");
		SPS1SMABox.SetWindowTextW(L"17512783.0");
		SPS2SMABox.SetWindowTextW(L"-390049424.0");
		TAtlanticBox.SetWindowTextW(L"1400.0");
		AtlanticLatitudeBox.SetWindowTextW(L"28.29028886");
		AtlanticLongitudeBox.SetWindowTextW(L"-19.5");
		TPacificBox.SetWindowTextW(L"30921.42");
		PacificLatitudeBox.SetWindowTextW(L"30.04649677");
		PacificLongitudeBox.SetWindowTextW(L"-171.0");

		//Temporary, until there is input boxes
		agc.BLOCKI.POLYCOFF[0] = 0.4788289;
		agc.BLOCKI.POLYCOFF[1] = 0.3411991e-1;
		agc.BLOCKI.POLYCOFF[2] = 0.1178593e-1;
		agc.BLOCKI.POLYCOFF[3] = 0.5318752e-4;
		agc.BLOCKI.POLYCOFF[4] = -0.2937101e-5;
		agc.BLOCKI.POLYCOFF[5] = 0.232939e-7;
		agc.BLOCKI.POLYCOFF[6] = -0.5793149e-10;
		break;
	}
}
