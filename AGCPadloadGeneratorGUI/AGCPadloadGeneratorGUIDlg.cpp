
// AGCPadloadGeneratorGUIDlg.cpp:
//

#include "pch.h"
#include "AGCPadloadGeneratorGUI.h"
#include "AGCPadloadGeneratorGUIDlg.h"
#include "afxdialogex.h"
#include <sstream>

#define CMC_COLOSSUS237 0
#define CMC_COLOSSUS249 1
#define CMC_COMANCE045 2
#define CMC_COMANCE055 3
#define CMC_COMANCE067 4
#define CMC_COMANCE072 5
#define CMC_COMANCE0108 6
#define CMC_ARTEMIS072 7
#define CMC_ARTEMIS072NBY70 8
#define CMC_ARTEMIS072NBY71 9

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
	DDX_Control(pDX, IDC_EDIT18, LATSPLBox);
	DDX_Control(pDX, IDC_EDIT19, LNGSPLBox);
	DDX_Control(pDX, IDC_EDIT20, CSMMASSBox);
	DDX_Control(pDX, IDC_EDIT21, LEMMASSBox);
	DDX_Control(pDX, IDC_EDIT22, PACTOFFBox);
	DDX_Control(pDX, IDC_EDIT52, YACTOFFBox);
	DDX_Control(pDX, IDC_EDIT23, LADPADBox);
	DDX_Control(pDX, IDC_EDIT24, LODPADBox);
	DDX_Control(pDX, IDC_EDIT25, ALFAPADBox);
	DDX_Control(pDX, IDC_EDIT26, P37RANGEBox);
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
	RopeNameBox.AddString(L"Comanche072");
	RopeNameBox.AddString(L"Comanche108");
	RopeNameBox.AddString(L"Artemis072");
	RopeNameBox.AddString(L"Artemis072NBY70");
	RopeNameBox.AddString(L"Artemis072NBY71");
	RopeNameBox.SetCurSel(CMC_COMANCE055);

	UpdateRopeSpecificEditFields();

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
	MissionBox.AddString(L"Apollo 16");
	MissionBox.AddString(L"Apollo 17");
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
	LATSPLBox.SetWindowText(L"26.48");
	LNGSPLBox.SetWindowText(L"-17.05");
	CSMMASSBox.SetWindowText(L"63386.7");
	LEMMASSBox.SetWindowText(L"33275.6");
	PACTOFFBox.SetWindowText(L"-1.541");
	YACTOFFBox.SetWindowText(L"1.321");
	LADPADBox.SetWindowText(L"0.3");
	LODPADBox.SetWindowText(L"0.18");
	ALFAPADBox.SetWindowText(L"-20.0");
	P37RANGEBox.SetWindowText(L"0.0");

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
	agc.BLOCKII.WRENDPOS = Utilities::Text2Double(&WRENDPOSBox);
	agc.BLOCKII.WRENDVEL = Utilities::Text2Double(&WRENDVELBox);
	agc.BLOCKII.RMAX = Utilities::Text2Double(&RMAXBox);
	agc.BLOCKII.VMAX = Utilities::Text2Double(&VMAXBox);
	agc.BLOCKII.LAT_SPL = Utilities::Text2Double(&LATSPLBox);
	agc.BLOCKII.LNG_SPL = Utilities::Text2Double(&LNGSPLBox);
	agc.BLOCKII.CSMMass = Utilities::Text2Double(&CSMMASSBox);
	agc.BLOCKII.LMMass = Utilities::Text2Double(&LEMMASSBox);
	agc.BLOCKII.PACTOFF = Utilities::Text2Double(&PACTOFFBox);
	agc.BLOCKII.YACTOFF = Utilities::Text2Double(&YACTOFFBox);
	agc.BLOCKII.LADPAD = Utilities::Text2Double(&LADPADBox);
	agc.BLOCKII.LODPAD = Utilities::Text2Double(&LODPADBox);
	agc.BLOCKII.ALFAPAD = Utilities::Text2Double(&ALFAPADBox);
	agc.BLOCKII.P37RANGE = Utilities::Text2Double(&P37RANGEBox);

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
	//UpdateData(FALSE);
}

void CAGCPadloadGeneratorGUIDlg::OnCbnSelchangeCombo2()
{
	UpdateRopeSpecificEditFields();
	//UpdateData(FALSE);
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
	LATSPLBox.SetWindowText(L"26.48");
	LNGSPLBox.SetWindowText(L"-17.05");

	switch (MissionBox.GetCurSel())
	{
	case 1: //Apollo 7
		RopeNameBox.SetCurSel(CMC_COLOSSUS237);
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
		CSMMASSBox.SetWindowTextW(L"32816.3");
		LEMMASSBox.SetWindowTextW(L"0.0");

		//TBD
		PACTOFFBox.SetWindowTextW(L"-1.541");
		YACTOFFBox.SetWindowTextW(L"1.321");
		LADPADBox.SetWindowTextW(L"0.27");
		LODPADBox.SetWindowTextW(L"0.207");
		ALFAPADBox.SetWindowTextW(L"-19.55");
		break;
	case 2: //Apollo 8
		RopeNameBox.SetCurSel(CMC_COLOSSUS237);
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
		CSMMASSBox.SetWindowTextW(L"64001.0");
		LEMMASSBox.SetWindowTextW(L"0.0");

		//TBD
		PACTOFFBox.SetWindowTextW(L"-1.541");
		YACTOFFBox.SetWindowTextW(L"1.321");
		LADPADBox.SetWindowTextW(L"0.27");
		LODPADBox.SetWindowTextW(L"0.207");
		ALFAPADBox.SetWindowTextW(L"-19.55");
		break;
	case 3: //Apollo 9
		RopeNameBox.SetCurSel(CMC_COLOSSUS249);
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
		CSMMASSBox.SetWindowTextW(L"59183.0");
		LEMMASSBox.SetWindowTextW(L"32000.0");

		//TBD
		PACTOFFBox.SetWindowTextW(L"-1.541");
		YACTOFFBox.SetWindowTextW(L"1.321");
		LADPADBox.SetWindowTextW(L"0.27");
		LODPADBox.SetWindowTextW(L"0.207");
		ALFAPADBox.SetWindowTextW(L"-19.55");
		break;
	case 4: //Apollo 10
		RopeNameBox.SetCurSel(CMC_COMANCE045);
		LaunchMJDInput.SetWindowTextW(L"40359.7006944");
		Launchpad.SetCurSel(1); //LC-39A
		RTEDBox.SetWindowTextW(L"1.6602637");
		LSAltitudeBox.SetWindowTextW(L"-3073.26");
		LSLatitudeBox.SetWindowTextW(L"0.71388");
		LSLongitudeBox.SetWindowTextW(L"23.707773");
		EMSAltBox.SetWindowTextW(L"294084.3");
		LaunchAzimuthBox.SetWindowTextW(L"72.028");
		HORIZALTBox.SetWindowTextW(L"24000");
		ALTVARBox.SetWindowText(L"1.5258e-05");
		CSMMASSBox.SetWindowTextW(L"63904.0");
		LEMMASSBox.SetWindowTextW(L"31579.7");
		PACTOFFBox.SetWindowTextW(L"-1.48");
		YACTOFFBox.SetWindowTextW(L"1.35");

		//TBD
		LADPADBox.SetWindowTextW(L"0.27");
		LODPADBox.SetWindowTextW(L"0.207");
		ALFAPADBox.SetWindowTextW(L"-19.55");
		P37RANGEBox.SetWindowTextW(L"1221.5");
		break;
	case 5: //Apollo 11
		RopeNameBox.SetCurSel(CMC_COMANCE055);
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
		CSMMASSBox.SetWindowTextW(L"63386.7");
		LEMMASSBox.SetWindowTextW(L"33275.6");
		PACTOFFBox.SetWindowTextW(L"-1.541");
		YACTOFFBox.SetWindowTextW(L"1.321");
		LADPADBox.SetWindowTextW(L"0.27");
		LODPADBox.SetWindowTextW(L"0.207");
		ALFAPADBox.SetWindowTextW(L"-19.55");
		P37RANGEBox.SetWindowTextW(L"1221.5");
		break;
	case 6: //Apollo 12
		Apollo12Padload();
		break;
	case 7: //Apollo 13
		Apollo13Padload();
		break;
	case 8: //Apollo 14
		Apollo14Padload();
		break;
	case 9: //Apollo 15
		Apollo15Padload();
		break;
	case 10: //Apollo 16
		Apollo16Padload();
		break;
	case 11: //Apollo 17
		Apollo17Padload();
		break;
	default:
		break;
	}

	UpdateRopeSpecificEditFields();
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

void CAGCPadloadGeneratorGUIDlg::UpdateRopeSpecificEditFields()
{
	if (RopeNameBox.GetCurSel() <= CMC_COLOSSUS249)
	{
		P37RANGEBox.SetReadOnly(true);
	}
	else
	{
		P37RANGEBox.SetReadOnly(false);
	}
}

void CAGCPadloadGeneratorGUIDlg::Apollo12Padload()
{
	RopeNameBox.SetCurSel(CMC_COMANCE067);
	EphemerisSpanBox.SetWindowTextW(L"14.5");
	LaunchMJDInput.SetWindowTextW(L"40539.6819444");
	Launchpad.SetCurSel(1); //LC-39A
	RTEDBox.SetWindowTextW(L"1.6602637");
	LSAltitudeBox.SetWindowTextW(L"-2371.27");
	LSLatitudeBox.SetWindowTextW(L"-2.9822165");
	LSLongitudeBox.SetWindowTextW(L"-23.391933");
	EMSAltBox.SetWindowTextW(L"304519.2");
	LaunchAzimuthBox.SetWindowTextW(L"72.029345");
	HORIZALTBox.SetWindowTextW(L"24000");
	CSMMASSBox.SetWindowTextW(L"63477.0");
	LEMMASSBox.SetWindowTextW(L"33559.3");
	PACTOFFBox.SetWindowTextW(L"-1.541");
	YACTOFFBox.SetWindowTextW(L"1.321");
	LADPADBox.SetWindowTextW(L"0.27");
	LODPADBox.SetWindowTextW(L"0.207");
	ALFAPADBox.SetWindowTextW(L"-20.5");
	P37RANGEBox.SetWindowTextW(L"1205.8");
}

void CAGCPadloadGeneratorGUIDlg::Apollo13Padload()
{
	RopeNameBox.SetCurSel(CMC_COMANCE055);
	EphemerisSpanBox.SetWindowTextW(L"14.5");
	LaunchMJDInput.SetWindowText(L"40687.8006944444444");
	Launchpad.SetCurSel(1); //LC-39A
	RTEDBox.SetWindowTextW(L"1.6602637");
	LSLatitudeBox.SetWindowText(L"-3.6686");
	LSLongitudeBox.SetWindowText(L"-17.4842");
	LSAltitudeBox.SetWindowText(L"-1405.0");
	EMSAltBox.SetWindowTextW(L"290001.0");
	LaunchAzimuthBox.SetWindowTextW(L"72.042916");
	HORIZALTBox.SetWindowTextW(L"24000");
	CSMMASSBox.SetWindowTextW(L"63678.1");
	LEMMASSBox.SetWindowTextW(L"33420.9");
	PACTOFFBox.SetWindowTextW(L"-1.517");
	YACTOFFBox.SetWindowTextW(L"1.3158");
	LADPADBox.SetWindowTextW(L"0.3");
	LODPADBox.SetWindowTextW(L"0.18");
	ALFAPADBox.SetWindowTextW(L"-21.49");
	P37RANGEBox.SetWindowTextW(L"1185.64");

	agc.BLOCKII.POLYNUM[0] = -7.646894e-2;
	agc.BLOCKII.POLYNUM[1] = 1.79988e-1;
	agc.BLOCKII.POLYNUM[2] = 9.298907e-4;
	agc.BLOCKII.POLYNUM[3] = -1.49242e-4;
	agc.BLOCKII.POLYNUM[4] = 2.078269e-6;
	agc.BLOCKII.POLYNUM[5] = -1.601873e-8;
	agc.BLOCKII.POLYNUM[6] = 4.401981e-11;
	agc.BLOCKII.RPSTART = 11.85;
	agc.BLOCKII.POLYSTOP = 149.5;
}

void CAGCPadloadGeneratorGUIDlg::Apollo14Padload()
{
	RopeNameBox.SetCurSel(CMC_ARTEMIS072NBY71);
	EphemerisSpanBox.SetWindowTextW(L"14.5");
	LaunchMJDInput.SetWindowTextW(L"40982.84930555555");
	Launchpad.SetCurSel(1); //LC-39A
	RTEDBox.SetWindowTextW(L"1.6602637");
	LSAltitudeBox.SetWindowTextW(L"-1405.2");
	LSLatitudeBox.SetWindowTextW(L"-3.67329493");
	LSLongitudeBox.SetWindowTextW(L"-17.46428902");
	EMSAltBox.SetWindowTextW(L"293597.2");
	LaunchAzimuthBox.SetWindowTextW(L"72.066946");
	HORIZALTBox.SetWindowTextW(L"28000");
	CSMMASSBox.SetWindowTextW(L"64457.1");
	LEMMASSBox.SetWindowTextW(L"33678.4");
	PACTOFFBox.SetWindowTextW(L"-1.416");
	YACTOFFBox.SetWindowTextW(L"1.314");
	LADPADBox.SetWindowTextW(L"0.27");
	LODPADBox.SetWindowTextW(L"0.207");
	ALFAPADBox.SetWindowTextW(L"-19.06");
	P37RANGEBox.SetWindowTextW(L"1187.35");

	agc.BLOCKII.POLYNUM[0] = 1.405176e-1;
	agc.BLOCKII.POLYNUM[1] = 2.2886283e-1;
	agc.BLOCKII.POLYNUM[2] = 4.4197154e-3;
	agc.BLOCKII.POLYNUM[3] = 1.099765e-5;
	agc.BLOCKII.POLYNUM[4] = -1.4904606e-7;
	agc.BLOCKII.POLYNUM[5] = -2.3591821e-9;
	agc.BLOCKII.POLYNUM[6] = 1.3334313e-11;
	agc.BLOCKII.RPSTART = 12.6;
	agc.BLOCKII.POLYSTOP = 150.0;
}

void CAGCPadloadGeneratorGUIDlg::Apollo15Padload()
{
	RopeNameBox.SetCurSel(CMC_ARTEMIS072);
	EphemerisSpanBox.SetWindowTextW(L"14.5");
	LaunchMJDInput.SetWindowTextW(L"41158.565277778");
	Launchpad.SetCurSel(1); //LC-39A
	RTEDBox.SetWindowTextW(L"1.6602637");
	LSAltitudeBox.SetWindowTextW(L"-3550.0");
	LSLatitudeBox.SetWindowTextW(L"26.07389");
	LSLongitudeBox.SetWindowTextW(L"3.65389");
	EMSAltBox.SetWindowTextW(L"300781.4");
	LaunchAzimuthBox.SetWindowTextW(L"80.08868");
	LATSPLBox.SetWindowText(L"20.3");
	LNGSPLBox.SetWindowText(L"-19.5");
	HORIZALTBox.SetWindowTextW(L"28000");
	CSMMASSBox.SetWindowTextW(L"66915.2");
	LEMMASSBox.SetWindowTextW(L"36206.2");
	PACTOFFBox.SetWindowTextW(L"-0.524");
	YACTOFFBox.SetWindowTextW(L"1.895");
	LADPADBox.SetWindowTextW(L"0.27");
	LODPADBox.SetWindowTextW(L"0.207");
	ALFAPADBox.SetWindowTextW(L"-18.64");
	P37RANGEBox.SetWindowTextW(L"1114.55");

	agc.BLOCKII.POLYNUM[0] = 1.7813387e-1;
	agc.BLOCKII.POLYNUM[1] = 1.8310605e-2;
	agc.BLOCKII.POLYNUM[2] = 1.1728888e-2;
	agc.BLOCKII.POLYNUM[3] = -1.1888579e-7;
	agc.BLOCKII.POLYNUM[4] = -1.642047e-6;
	agc.BLOCKII.POLYNUM[5] = 1.270428e-8;
	agc.BLOCKII.POLYNUM[6] = -2.9239175e-11;
	agc.BLOCKII.RPSTART = 10.97;
	agc.BLOCKII.POLYSTOP = 144.25;
}

void CAGCPadloadGeneratorGUIDlg::Apollo16Padload()
{
	RopeNameBox.SetCurSel(CMC_ARTEMIS072);
	EphemerisSpanBox.SetWindowTextW(L"14.5");
	LaunchMJDInput.SetWindowTextW(L"41423.7458333333");
	Launchpad.SetCurSel(1); //LC-39A
	RTEDBox.SetWindowTextW(L"1.6602637");
	LSAltitudeBox.SetWindowTextW(L"-260.0");
	LSLatitudeBox.SetWindowTextW(L"-9.00028");
	LSLongitudeBox.SetWindowTextW(L"15.516389");
	EMSAltBox.SetWindowTextW(L"290000.1");
	LaunchAzimuthBox.SetWindowTextW(L"72.03443909");
	HORIZALTBox.SetWindowTextW(L"28000");
	CSMMASSBox.SetWindowTextW(L"66887.6");
	LEMMASSBox.SetWindowTextW(L"36208.3");
	PACTOFFBox.SetWindowTextW(L"-0.527");
	YACTOFFBox.SetWindowTextW(L"1.898");
	LADPADBox.SetWindowTextW(L"0.27");
	LODPADBox.SetWindowTextW(L"0.207");
	ALFAPADBox.SetWindowTextW(L"-18.51");
	P37RANGEBox.SetWindowTextW(L"1086.533");
	LATSPLBox.SetWindowText(L"26.5");
	LNGSPLBox.SetWindowText(L"-17.0");

	agc.BLOCKII.POLYNUM[0] = -9.5207497e-2;
	agc.BLOCKII.POLYNUM[1] = 1.9666296e-1;
	agc.BLOCKII.POLYNUM[2] = 9.5150055e-3;
	agc.BLOCKII.POLYNUM[3] = -1.4498524e-4;
	agc.BLOCKII.POLYNUM[4] = 1.9137262e-6;
	agc.BLOCKII.POLYNUM[5] = -1.5121426e-8;
	agc.BLOCKII.POLYNUM[6] = 4.3469835e-11;
	agc.BLOCKII.RPSTART = 11.85;
	agc.BLOCKII.POLYSTOP = 147.0;
}

void CAGCPadloadGeneratorGUIDlg::Apollo17Padload()
{
	RopeNameBox.SetCurSel(CMC_ARTEMIS072);
	EphemerisSpanBox.SetWindowTextW(L"14.5");
	LaunchMJDInput.SetWindowTextW(L"41658.120138888");
	Launchpad.SetCurSel(1); //LC-39A
	RTEDBox.SetWindowTextW(L"1.6602637");
	LSAltitudeBox.SetWindowTextW(L"-3606.0");
	LSLatitudeBox.SetWindowTextW(L"20.164029");
	LSLongitudeBox.SetWindowTextW(L"30.749532");
	EMSAltBox.SetWindowTextW(L"290000.1");
	LaunchAzimuthBox.SetWindowTextW(L"72.141393");
	HORIZALTBox.SetWindowTextW(L"28000");
	CSMMASSBox.SetWindowTextW(L"66893.4");
	LEMMASSBox.SetWindowTextW(L"36255.4");
	PACTOFFBox.SetWindowTextW(L"-0.578");
	YACTOFFBox.SetWindowTextW(L"1.892");
	LADPADBox.SetWindowTextW(L"0.27");
	LODPADBox.SetWindowTextW(L"0.207");
	ALFAPADBox.SetWindowTextW(L"-18.97");
	P37RANGEBox.SetWindowTextW(L"1074.63");
	LATSPLBox.SetWindowText(L"26.35");
	LNGSPLBox.SetWindowText(L"-17.1");

	agc.BLOCKII.POLYNUM[0] = -3.7352617e-1;
	agc.BLOCKII.POLYNUM[1] = 2.4436145e-1;
	agc.BLOCKII.POLYNUM[2] = 6.2373017e-3;
	agc.BLOCKII.POLYNUM[3] = -9.0587464e-5;
	agc.BLOCKII.POLYNUM[4] = 1.685017e-6;
	agc.BLOCKII.POLYNUM[5] = -1.6095034e-8;
	agc.BLOCKII.POLYNUM[6] = 4.9817897e-11;
	agc.BLOCKII.RPSTART = 12.25;
	agc.BLOCKII.POLYSTOP = 147.75;
}