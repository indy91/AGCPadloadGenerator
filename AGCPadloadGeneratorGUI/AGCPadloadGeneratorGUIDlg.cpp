
// AGCPadloadGeneratorGUIDlg.cpp:
//

#include "pch.h"
#include "AGCPadloadGeneratorGUI.h"
#include "AGCPadloadGeneratorGUIDlg.h"
#include "afxdialogex.h"
#include <sstream>

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

// CAGCPadloadGeneratorGUIDlg


CAGCPadloadGeneratorGUIDlg::CAGCPadloadGeneratorGUIDlg(CWnd* pParent /*=nullptr*/)
	: CDialogEx(IDD_AGCPADLOADGENERATORGUI_DIALOG, pParent)
	, LaunchPadValue(_T(""))
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CAGCPadloadGeneratorGUIDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_EDIT1, LaunchMJDInput);
	DDX_Control(pDX, IDC_COMBO1, Launchpad);
	DDX_CBString(pDX, IDC_COMBO1, LaunchPadValue);
	DDX_Control(pDX, IDC_EDIT2, TLANDBox);
	DDX_Control(pDX, IDC_EDIT3, LSLatitudeBox);
	DDX_Control(pDX, IDC_EDIT4, LSLongitudeBox);
	DDX_Control(pDX, IDC_EDIT5, LSAltitudeBox);
	DDX_Control(pDX, IDC_EDIT6, RTEDBox);
	DDX_Control(pDX, IDC_COMBO2, RopeNameBox);
	DDX_Control(pDX, IDC_EDIT7, EMSAltBox);
	DDX_Control(pDX, IDC_EDIT8, UNITWBox);
	DDX_Control(pDX, IDC_EDIT9, EphemerisSpanBox);
	DDX_Control(pDX, IDC_COMBO3, MissionBox);
	DDX_Control(pDX, IDC_EDIT10, LaunchAzimuthBox);
	DDX_Control(pDX, IDC_EDIT11, CDUCHKWDBox);
	DDX_Control(pDX, IDC_EDIT12, HORIZALTBox);
	DDX_Control(pDX, IDC_EDIT13, ALTVARBox);
	DDX_Control(pDX, IDC_EDIT14, WRENDPOSBox);
	DDX_Control(pDX, IDC_EDIT15, WRENDVELBox);
	DDX_Control(pDX, IDC_EDIT16, RMAXBox);
	DDX_Control(pDX, IDC_EDIT17, VMAXBox);
}

BEGIN_MESSAGE_MAP(CAGCPadloadGeneratorGUIDlg, CDialogEx)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	ON_BN_CLICKED(IDOK, &CAGCPadloadGeneratorGUIDlg::OnBnClickedOk)
	ON_CBN_SELCHANGE(IDC_COMBO1, &CAGCPadloadGeneratorGUIDlg::OnCbnSelchangeCombo1)
	ON_CBN_SELCHANGE(IDC_COMBO2, &CAGCPadloadGeneratorGUIDlg::OnCbnSelchangeCombo2)
	ON_CBN_SELCHANGE(IDC_COMBO3, &CAGCPadloadGeneratorGUIDlg::OnCbnSelchangeCombo3)
END_MESSAGE_MAP()

BOOL CAGCPadloadGeneratorGUIDlg::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);
	ASSERT(IDM_ABOUTBOX < 0xF000);

	CMenu* pSysMenu = GetSystemMenu(FALSE);
	if (pSysMenu != nullptr)
	{
		BOOL bNameValid;
		CString strAboutMenu;
		bNameValid = strAboutMenu.LoadString(IDS_ABOUTBOX);
		ASSERT(bNameValid);
		if (!strAboutMenu.IsEmpty())
		{
			pSysMenu->AppendMenu(MF_SEPARATOR);
			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);
		}
	}

	SetIcon(m_hIcon, TRUE);
	SetIcon(m_hIcon, FALSE);

	Launchpad.AddString(L"LC-34");
	Launchpad.AddString(L"LC-39A");
	Launchpad.AddString(L"LC-39B");
	Launchpad.SetCurSel(1);

	RopeNameBox.AddString(L"Colossus237");
	RopeNameBox.AddString(L"Colossus249");
	RopeNameBox.AddString(L"Comanche045");
	RopeNameBox.AddString(L"Comanche055");
	RopeNameBox.AddString(L"Comanche067");
	RopeNameBox.AddString(L"Artemis072");
	RopeNameBox.SetCurSel(3);

	MissionBox.AddString(L"Manual");
	MissionBox.AddString(L"Apollo 7");
	MissionBox.AddString(L"Apollo 8");
	MissionBox.AddString(L"Apollo 9");
	MissionBox.AddString(L"Apollo 10");
	MissionBox.AddString(L"Apollo 11");
	MissionBox.AddString(L"Apollo 12");
	MissionBox.AddString(L"Apollo 13");
	MissionBox.AddString(L"Apollo 14");
	MissionBox.AddString(L"Apollo 15");
	MissionBox.SetCurSel(0);

	TLANDBox.SetWindowText(L"4.5");
	UNITWBox.SetWindowText(L"-0.5");
	LSLatitudeBox.SetWindowText(L"0.0");
	LSLongitudeBox.SetWindowText(L"0.0");
	LSAltitudeBox.SetWindowText(L"0.0");
	RTEDBox.SetWindowText(L"1.69107");//1.6602637");
	LaunchMJDInput.SetWindowText(L"0.0");
	EMSAltBox.SetWindowText(L"284643");//290000.0");
	EphemerisSpanBox.SetWindowText(L"14.5");
	LaunchAzimuthBox.SetWindowText(L"72.0");
	CDUCHKWDBox.SetWindowText(L"5");
	HORIZALTBox.SetWindowText(L"28000");
	ALTVARBox.SetWindowText(L"1.52168e-5");
	WRENDPOSBox.SetWindowText(L"3048");
	WRENDVELBox.SetWindowText(L"3.048");
	RMAXBox.SetWindowText(L"2000");
	VMAXBox.SetWindowText(L"2");

	return TRUE;
}

void CAGCPadloadGeneratorGUIDlg::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
	}
	else
	{
		CDialogEx::OnSysCommand(nID, lParam);
	}
}

void CAGCPadloadGeneratorGUIDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this);

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialogEx::OnPaint();
	}
}

HCURSOR CAGCPadloadGeneratorGUIDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}

void CAGCPadloadGeneratorGUIDlg::OnBnClickedOk()
{
	CString string;

	agc.LaunchMJD = Utilities::Text2Double(&LaunchMJDInput);
	agc.T_504LM = Utilities::Text2Double(&TLANDBox);
	agc.T_UNITW = Utilities::Text2Double(&UNITWBox);
	agc.LSLat = Utilities::Text2Double(&LSLatitudeBox);
	agc.LSLng = Utilities::Text2Double(&LSLongitudeBox);
	agc.LSAlt = Utilities::Text2Double(&LSAltitudeBox);
	agc.RTED1 = Utilities::Text2Double(&RTEDBox);
	agc.EMSALT = Utilities::Text2Double(&EMSAltBox);
	agc.EphemerisSpan = Utilities::Text2Double(&EphemerisSpanBox);
	agc.LaunchAzimuth = Utilities::Text2Double(&LaunchAzimuthBox);
	agc.CDUCHKWD = Utilities::Text2Double(&CDUCHKWDBox) / 100.0;
	agc.HORIZALT = Utilities::Text2Double(&HORIZALTBox);
	agc.ALTVAR = Utilities::Text2Double(&ALTVARBox);
	agc.WRENDPOS = Utilities::Text2Double(&WRENDPOSBox);
	agc.WRENDVEL = Utilities::Text2Double(&WRENDVELBox);
	agc.RMAX = Utilities::Text2Double(&RMAXBox)*0.3048;
	agc.VMAX = Utilities::Text2Double(&VMAXBox)*0.3048;

	Launchpad.GetWindowText(string);
	std::wstring ws = std::wstring(string.GetString());
	agc.Pad = std::string(ws.begin(), ws.end());

	RopeNameBox.GetWindowText(string);
	ws = std::wstring(string.GetString());
	agc.RopeName = std::string(ws.begin(), ws.end());

	agc.RunCMC();
}

void CAGCPadloadGeneratorGUIDlg::OnCbnSelchangeCombo1()
{
	Launchpad.GetLBText(Launchpad.GetCurSel(), LaunchPadValue);
	UpdateData(FALSE);
}

void CAGCPadloadGeneratorGUIDlg::OnCbnSelchangeCombo2()
{
	RopeNameBox.GetLBText(RopeNameBox.GetCurSel(), RopeNameValue);
	UpdateData(FALSE);
}

void CAGCPadloadGeneratorGUIDlg::OnCbnSelchangeCombo3()
{
	MissionBox.GetLBText(MissionBox.GetCurSel(), MissionNameValue);
	UpdateData(FALSE);

	agc.RPVAR = 2000.0*2000.0;
	agc.S22WSUBL = 10000.0;

	CDUCHKWDBox.SetWindowTextW(L"5");
	ALTVARBox.SetWindowText(L"1.52168e-5");
	WRENDPOSBox.SetWindowText(L"3048");
	WRENDVELBox.SetWindowText(L"3.048");
	RMAXBox.SetWindowText(L"2000");
	VMAXBox.SetWindowText(L"2");

	switch (MissionBox.GetCurSel())
	{
	case 1: //Apollo 7
		RopeNameBox.SetCurSel(0); //Colossus 237
		EphemerisSpanBox.SetWindowTextW(L"14.5");
		LaunchMJDInput.SetWindowTextW(L"40140.62690972");
		Launchpad.SetCurSel(0); //LC-34
		RTEDBox.SetWindowTextW(L"1.69107");
		LSAltitudeBox.SetWindowTextW(L"0.0");
		LSLatitudeBox.SetWindowTextW(L"0.0");
		LSLongitudeBox.SetWindowTextW(L"0.0");
		EMSAltBox.SetWindowTextW(L"284643");
		LaunchAzimuthBox.SetWindowTextW(L"72.0");
		agc.S22WSUBL = 10000.0*0.3048;
		agc.RPVAR = pow(1500.0*0.3048, 2);
		CDUCHKWDBox.SetWindowTextW(L"0");
		HORIZALTBox.SetWindowTextW(L"18000");
		ALTVARBox.SetWindowText(L"1e-6");
		WRENDPOSBox.SetWindowText(L"1000");
		WRENDVELBox.SetWindowText(L"1");
		RMAXBox.SetWindowText(L"-1");
		VMAXBox.SetWindowText(L"-1");
		break;
	case 2: //Apollo 8
		RopeNameBox.SetCurSel(0); //Colossus 237
		EphemerisSpanBox.SetWindowTextW(L"10.5");
		LaunchMJDInput.SetWindowTextW(L"40211.53541666666");
		Launchpad.SetCurSel(1); //LC-39A
		RTEDBox.SetWindowTextW(L"1.6602637");
		LSAltitudeBox.SetWindowTextW(L"-1518.64");
		LSLatitudeBox.SetWindowTextW(L"2.6317");
		LSLongitudeBox.SetWindowTextW(L"34.0253");
		EMSAltBox.SetWindowTextW(L"297431");
		LaunchAzimuthBox.SetWindowTextW(L"72.124");
		HORIZALTBox.SetWindowTextW(L"18000");
		ALTVARBox.SetWindowText(L"1e-6");
		WRENDPOSBox.SetWindowText(L"1000");
		WRENDVELBox.SetWindowText(L"1");
		break;
	case 3: //Apollo 9
		RopeNameBox.SetCurSel(1); //Colossus 249
		LaunchMJDInput.SetWindowTextW(L"40283.666667");
		Launchpad.SetCurSel(1); //LC-39A
		RTEDBox.SetWindowTextW(L"1.69107");
		EMSAltBox.SetWindowTextW(L"284643");
		LaunchAzimuthBox.SetWindowTextW(L"72.0");
		CDUCHKWDBox.SetWindowTextW(L"0");
		HORIZALTBox.SetWindowTextW(L"18000");
		ALTVARBox.SetWindowText(L"1e-6");
		WRENDPOSBox.SetWindowText(L"1000");
		WRENDVELBox.SetWindowText(L"1");
		break;
	case 4: //Apollo 10
		HORIZALTBox.SetWindowTextW(L"24000");
		ALTVARBox.SetWindowText(L"1.5258e-05");
		break;
	case 5: //Apollo 11
		RopeNameBox.SetCurSel(3); //Comanche 55
		EphemerisSpanBox.SetWindowTextW(L"10.5");
		LaunchMJDInput.SetWindowTextW(L"40418.5638889");
		Launchpad.SetCurSel(1); //LC-39A
		RTEDBox.SetWindowTextW(L"1.6602637");
		LSAltitudeBox.SetWindowTextW(L"-2671.9684");
		LSLatitudeBox.SetWindowTextW(L"0.691395");
		LSLongitudeBox.SetWindowTextW(L"23.71689");
		EMSAltBox.SetWindowTextW(L"294084.3");
		LaunchAzimuthBox.SetWindowTextW(L"72.05897");
		HORIZALTBox.SetWindowTextW(L"24000");
		ALTVARBox.SetWindowText(L"1.5258e-05");
		break;
	case 6: //Apollo 12
		RopeNameBox.SetCurSel(4); //Comanche 67
		EphemerisSpanBox.SetWindowTextW(L"14.5");
		LaunchMJDInput.SetWindowTextW(L"40539.6819444");
		Launchpad.SetCurSel(1); //LC-39A
		RTEDBox.SetWindowTextW(L"1.6602637");
		LSAltitudeBox.SetWindowTextW(L"-2371.27");
		LSLatitudeBox.SetWindowTextW(L"-2.9822165");
		LSLongitudeBox.SetWindowTextW(L"-23.391933");
		EMSAltBox.SetWindowTextW(L"297431.0");
		LaunchAzimuthBox.SetWindowTextW(L"72.029345");
		HORIZALTBox.SetWindowTextW(L"24000");
		break;
	case 7: //Apollo 13
		Launchpad.SetCurSel(1); //LC-39A
		RTEDBox.SetWindowTextW(L"1.6602637");
		HORIZALTBox.SetWindowTextW(L"24000");
		break;
	case 8: //Apollo 14
		Launchpad.SetCurSel(1); //LC-39A
		RTEDBox.SetWindowTextW(L"1.6602637");
		HORIZALTBox.SetWindowTextW(L"28000");
		break;
	case 9: //Apollo 15
		RopeNameBox.SetCurSel(5); //Artemis 72
		EphemerisSpanBox.SetWindowTextW(L"14.5");
		LaunchMJDInput.SetWindowTextW(L"41158.565277778");
		Launchpad.SetCurSel(1); //LC-39A
		RTEDBox.SetWindowTextW(L"1.6602637");
		LSAltitudeBox.SetWindowTextW(L"-3550.284");
		LSLatitudeBox.SetWindowTextW(L"26.074");
		LSLongitudeBox.SetWindowTextW(L"3.654");
		EMSAltBox.SetWindowTextW(L"297431.0");
		LaunchAzimuthBox.SetWindowTextW(L"80.08868");
		HORIZALTBox.SetWindowTextW(L"28000");
		break;
	case 10: //Apollo 16
		HORIZALTBox.SetWindowTextW(L"28000");
		break;
	case 11: //Apollo 17
		HORIZALTBox.SetWindowTextW(L"28000");
		break;
	default:
		break;
	}
}

void CAGCPadloadGeneratorGUIDlg::Double2Text(double val, CEdit *ed, int length)
{
	CString string;
	string.Format(_T("%.*lf"), length, val);
	ed->SetWindowText(string);
}

void CAGCPadloadGeneratorGUIDlg::String2Text(std::string val, CEdit *ed)
{
	CString string(val.c_str());
	ed->SetWindowText(string);
}

void CAGCPadloadGeneratorGUIDlg::String2Text(std::string val, CComboBox *ed)
{
	CString string(val.c_str());
	ed->SetWindowText(string);
}