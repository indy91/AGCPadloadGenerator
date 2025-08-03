
// AGCPadloadGeneratorGUIDlg.cpp:
//

#include "pch.h"
#include "AGCPadloadGeneratorGUI.h"
#include "AGCPadloadGeneratorGUIDlg.h"
#include "afxdialogex.h"
#include <sstream>

#define CMC_COLOSSUS237 0
#define CMC_COLOSSUS249 1
#define CMC_COMANCHE045 2
#define CMC_COMANCHE055 3
#define CMC_COMANCHE067 4
#define CMC_COMANCHE072 5
#define CMC_COMANCHE108 6
#define CMC_ARTEMIS072 7

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

// CAGCPadloadGeneratorGUIDlg


CAGCPadloadGeneratorGUIDlg::CAGCPadloadGeneratorGUIDlg(CWnd* pParent /*=nullptr*/)
	: CDialogEx(IDD_AGCPADLOADGENERATORGUI_DIALOG, pParent)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CAGCPadloadGeneratorGUIDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_EDIT1, LaunchMJDInput);
	DDX_Control(pDX, IDC_COMBO1, Launchpad);
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
	DDX_Control(pDX, IDC_EDIT27, EMDOTBox);
	DDX_Control(pDX, IDC_EDIT28, MinImp1Box);
	DDX_Control(pDX, IDC_EDIT29, MinImp2Box);
	DDX_Control(pDX, IDC_EDIT30, MinImp3Box);
	DDX_Control(pDX, IDC_EDIT31, MinImp4Box);
	DDX_Control(pDX, IDC_MINIMP_1, MinImp1Label);
	DDX_Control(pDX, IDC_MINIMP_2, MinImp2Label);
	DDX_Control(pDX, IDC_MINIMP_3, MinImp3Label);
	DDX_Control(pDX, IDC_MINIMP_4, MinImp4Label);
	DDX_Control(pDX, IDC_MINIMP_1U, MinImp1Unit);
	DDX_Control(pDX, IDC_MINIMP_2U, MinImp2Unit);
	DDX_Control(pDX, IDC_MINIMP_3U, MinImp3Unit);
	DDX_Control(pDX, IDC_MINIMP_4U, MinImp4Unit);
	DDX_Control(pDX, IDC_EDIT32, TRUNSFBox);
	DDX_Control(pDX, IDC_EDIT33, SHAFTSFBox);
	DDX_Control(pDX, IDC_EDIT34, WMIDPOSBox);
	DDX_Control(pDX, IDC_EDIT35, WMIDVELBox);
	DDX_Control(pDX, IDC_EDIT36, RVARMINBox);
	DDX_Control(pDX, IDC_EDIT37, OutputBox);
	DDX_Control(pDX, IDC_EDIT38, PIOSDataSetBox);
	DDX_Control(pDX, IDC_CHECK2, R2ModelBox);
	DDX_Control(pDX, IDC_EDIT39, PBIASXBox);
	DDX_Control(pDX, IDC_EDIT40, PIPASCFXBox);
	DDX_Control(pDX, IDC_EDIT41, PBIASYBox);
	DDX_Control(pDX, IDC_EDIT42, PIPASCFYBox);
	DDX_Control(pDX, IDC_EDIT43, PBIASZBox);
	DDX_Control(pDX, IDC_EDIT44, PIPASCFZBox);
	DDX_Control(pDX, IDC_EDIT45, NBDXBox);
	DDX_Control(pDX, IDC_EDIT46, NBDYBox);
	DDX_Control(pDX, IDC_EDIT47, NBDZBox);
	DDX_Control(pDX, IDC_EDIT48, ADIAXBox);
	DDX_Control(pDX, IDC_EDIT49, ADIAYBox);
	DDX_Control(pDX, IDC_EDIT50, ADIAZBox);
	DDX_Control(pDX, IDC_EDIT51, ADSRAXBox);
	DDX_Control(pDX, IDC_EDIT65, ADSRAYBox);
	DDX_Control(pDX, IDC_EDIT66, ADSRAZBox);
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

	//Create the ToolTip control
	if (!m_ToolTip.Create(this))
	{
		TRACE0("Unable to create the ToolTip!");
	}
	else
	{
		// Add tool tips to the controls, either by hard coded string 
		// or using the string table resource

		m_ToolTip.AddTool(&PIOSDataSetBox, _T("Planetary Inertial Orientation Subroutine data from PIOSDataSets.txt"));

		m_ToolTip.AddTool(&UNITWBox, _T("Time in days from liftoff at which Earth rotations are exactly accurate"));
		m_ToolTip.AddTool(&TLANDBox, _T("Time in days from liftoff at which Moon rotations are exactly accurate"));
		m_ToolTip.AddTool(&EphemerisSpanBox, _T("Duration of accuracy of lunar ephemeris"));
		m_ToolTip.AddTool(&CDUCHKWDBox, _T("Delay before checking IMU CDU after an optics mark"));
		m_ToolTip.AddTool(&HORIZALTBox, _T("Horizon altitude used in P23 above Fischer ellipse"));
		m_ToolTip.AddTool(&RTEDBox, _T("Constant term for desired flight path angle cotangent in P37"));
		m_ToolTip.AddTool(&EMSAltBox, _T("Entry monitoring system initialization altitude"));
		m_ToolTip.AddTool(&LADPADBox, _T("Nominal L/D"));
		m_ToolTip.AddTool(&LODPADBox, _T("Final phase L/D"));
		m_ToolTip.AddTool(&ALFAPADBox, _T("Nominal entry trim angle (expected to be a negative number)"));
		m_ToolTip.AddTool(&P37RANGEBox, _T("Override on entry range (not time) in P37 if cell is non-zero"));
		m_ToolTip.AddTool(&PACTOFFBox, _T("TVC DAP pitch trim"));
		m_ToolTip.AddTool(&YACTOFFBox, _T("TVC DAP yaw trim"));
		m_ToolTip.AddTool(&EMDOTBox, _T("Mass flow rate for SPS engine"));
		m_ToolTip.AddTool(&LATSPLBox, _T("Target latitude for Mode III abort"));
		m_ToolTip.AddTool(&LNGSPLBox, _T("Target longitude for Mode III abort"));

		m_ToolTip.AddTool(&ALTVARBox, _T("A priori measurement accuracy of back-up optics in R23"));
		m_ToolTip.AddTool(&WRENDPOSBox, _T("W matrix position initialization for rendezvous"));
		m_ToolTip.AddTool(&WRENDVELBox, _T("W matrix velocity initialization for rendezvous"));
		m_ToolTip.AddTool(&RMAXBox, _T("Maximum automatic rendezvous position update"));
		m_ToolTip.AddTool(&VMAXBox, _T("Maximum automatic rendezvous velocity update"));
		m_ToolTip.AddTool(&RVARMINBox, _T("Minimum variance for VHF measurement"));
		m_ToolTip.AddTool(&WMIDPOSBox, _T("W-matrix position error initialization for midcourse navigation (P23)"));
		m_ToolTip.AddTool(&WMIDVELBox, _T("W-matrix velocity error initialization for midcourse navigation (P23)"));

		m_ToolTip.AddTool(&MinImp1Box, _T("None"));
		m_ToolTip.AddTool(&MinImp2Box, _T("None"));
		m_ToolTip.AddTool(&MinImp3Box, _T("None"));
		m_ToolTip.AddTool(&MinImp4Box, _T("None"));

		m_ToolTip.AddTool(&TRUNSFBox, _T("Scale factor for rate command to optics trunnion in P24"));
		m_ToolTip.AddTool(&SHAFTSFBox, _T("Scale factor for rate command to optics shaft in P24"));

		m_ToolTip.AddTool(&R2ModelBox, _T("R2 gravity model only supported by Open Orbiter"));

		m_ToolTip.Activate(TRUE);
	}

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
	RopeNameBox.SetCurSel(CMC_COMANCHE055);

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
	MissionBox.AddString(L"Apollo 11 July 18");
	MissionBox.AddString(L"Apollo 11 July 21");
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
	RVARMINBox.SetWindowText(L"40000");
	WMIDPOSBox.SetWindowText(L"30000");
	WMIDVELBox.SetWindowText(L"30");
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
	EMDOTBox.SetWindowText(L"65.272");

	TRUNSFBox.SetWindowText(L"0.0");
	SHAFTSFBox.SetWindowText(L"0.0");
	MinImp1Box.SetWindowText(L"0.0");
	MinImp2Box.SetWindowText(L"0.0");
	MinImp3Box.SetWindowText(L"0.0");
	MinImp4Box.SetWindowText(L"0.0");

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

BOOL CAGCPadloadGeneratorGUIDlg::PreTranslateMessage(MSG* pMsg)
{
	m_ToolTip.RelayEvent(pMsg);

	return CDialog::PreTranslateMessage(pMsg);
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
	agc.BLOCKII.RVARMIN = Utilities::Text2Double(&RVARMINBox);
	agc.CMCDATA.WMIDPOS = Utilities::Text2Double(&WMIDPOSBox);
	agc.CMCDATA.WMIDVEL = Utilities::Text2Double(&WMIDVELBox);
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
	agc.CMCDATA.EMDOT = Utilities::Text2Double(&EMDOTBox);

	agc.CMCDATA.TRUNSF = Utilities::Text2Double(&TRUNSFBox);
	agc.CMCDATA.SHAFTSF = Utilities::Text2Double(&SHAFTSFBox);

	if (RopeNameBox.GetCurSel() <= CMC_COMANCHE108)
	{
		agc.CMCDATA.EK1VAL = Utilities::Text2Double(&MinImp1Box);
		agc.CMCDATA.EK2VAL = Utilities::Text2Double(&MinImp2Box);
		agc.CMCDATA.EK3VAL = Utilities::Text2Double(&MinImp3Box);
		agc.CMCDATA.FANG = Utilities::Text2Double(&MinImp4Box);
	}
	else
	{
		agc.CMCDATA.EIMP1SEC = Utilities::Text2Double(&MinImp1Box);
		agc.CMCDATA.EFIMP01 = Utilities::Text2Double(&MinImp2Box);
		agc.CMCDATA.EFIMP16 = Utilities::Text2Double(&MinImp3Box);
	}

	Launchpad.GetWindowText(string);
	std::wstring ws = std::wstring(string.GetString());
	agc.Pad = std::string(ws.begin(), ws.end());

	RopeNameBox.GetWindowText(string);
	ws = std::wstring(string.GetString());
	agc.RopeName = std::string(ws.begin(), ws.end());

	PIOSDataSetBox.GetWindowText(string);
	ws = std::wstring(string.GetString());
	agc.PIOSDataSetName = std::string(ws.begin(), ws.end());

	agc.BLOCKII.R2Model = (R2ModelBox.GetCheck() != 0);

	agc.CMCDATA.IMUBiasCompensation.PBIASX = Utilities::Text2Double(&PBIASXBox);
	agc.CMCDATA.IMUBiasCompensation.PIPASCFX = Utilities::Text2Double(&PIPASCFXBox);
	agc.CMCDATA.IMUBiasCompensation.PBIASY = Utilities::Text2Double(&PBIASYBox);
	agc.CMCDATA.IMUBiasCompensation.PIPASCFY = Utilities::Text2Double(&PIPASCFYBox);
	agc.CMCDATA.IMUBiasCompensation.PBIASZ = Utilities::Text2Double(&PBIASZBox);
	agc.CMCDATA.IMUBiasCompensation.PIPASCFZ = Utilities::Text2Double(&PIPASCFZBox);
	agc.CMCDATA.IMUBiasCompensation.NBDX = Utilities::Text2Double(&NBDXBox);
	agc.CMCDATA.IMUBiasCompensation.NBDY = Utilities::Text2Double(&NBDYBox);
	agc.CMCDATA.IMUBiasCompensation.NBDZ = Utilities::Text2Double(&NBDZBox);
	agc.CMCDATA.IMUBiasCompensation.ADIAX = Utilities::Text2Double(&ADIAXBox);
	agc.CMCDATA.IMUBiasCompensation.ADIAY = Utilities::Text2Double(&ADIAYBox);
	agc.CMCDATA.IMUBiasCompensation.ADIAZ = Utilities::Text2Double(&ADIAZBox);
	agc.CMCDATA.IMUBiasCompensation.ADSRAX = Utilities::Text2Double(&ADSRAXBox);
	agc.CMCDATA.IMUBiasCompensation.ADSRAY = Utilities::Text2Double(&ADSRAYBox);
	agc.CMCDATA.IMUBiasCompensation.ADSRAZ = Utilities::Text2Double(&ADSRAZBox);

	int message = agc.RunCMC();

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

void CAGCPadloadGeneratorGUIDlg::OnCbnSelchangeCombo1()
{
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

	//Load some defaults
	agc.RPVAR = 2000.0*2000.0;
	agc.S22WSUBL = 10000.0;

	CDUCHKWDBox.SetWindowTextW(L"5");
	ALTVARBox.SetWindowText(L"1.52168e-5");
	WRENDPOSBox.SetWindowText(L"3048");
	WRENDVELBox.SetWindowText(L"3.048");
	RMAXBox.SetWindowText(L"2000");
	VMAXBox.SetWindowText(L"2");
	RVARMINBox.SetWindowText(L"40000");
	WMIDPOSBox.SetWindowText(L"30000");
	WMIDVELBox.SetWindowText(L"30");

	LATSPLBox.SetWindowText(L"26.48");
	LNGSPLBox.SetWindowText(L"-17.05");

	//Select new rope first
	switch (MissionBox.GetCurSel())
	{
	case 1: //Apollo 7
	case 2: //Apollo 8
		RopeNameBox.SetCurSel(CMC_COLOSSUS237);
		break;
	case 3: //Apollo 9
		RopeNameBox.SetCurSel(CMC_COLOSSUS249);
		break;
	case 4: //Apollo 10
		RopeNameBox.SetCurSel(CMC_COMANCHE045);
		break;
	case 5: //Apollo 11
		RopeNameBox.SetCurSel(CMC_COMANCHE055);
		break;
	case 6: //Apollo 12
		RopeNameBox.SetCurSel(CMC_COMANCHE067);
		break;
	case 7: //Apollo 13
		RopeNameBox.SetCurSel(CMC_COMANCHE072);
		break;
	case 8: //Apollo 14
		RopeNameBox.SetCurSel(CMC_ARTEMIS072); //TBD: Comanche 108
		break;
	case 9: //Apollo 15
	case 10: //Apollo 16
	case 11: //Apollo 17
		RopeNameBox.SetCurSel(CMC_ARTEMIS072);
		break;
	}

	//Update to new rope here, so mission specific numbers can be overwritten later
	UpdateRopeSpecificEditFields();

	//Now mission specific numbers
	switch (MissionBox.GetCurSel())
	{
	case 1: //Apollo 7
		Apollo7Padload();
		break;
	case 2: //Apollo 8
		Apollo8Padload();
		break;
	case 3: //Apollo 9
		Apollo9Padload();
		break;
	case 4: //Apollo 10
		Apollo10Padload();
		break;
	case 5: //Apollo 11
		Apollo11Padload(16);
		break;
	case 12: //Apollo 11 July 18 launch
		Apollo11Padload(18);
		break;
	case 13: //Apollo 11 July 21 launch
		Apollo11Padload(21);
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

	//Special cases for missions that currently use a rope created for a different launch year
	switch (MissionBox.GetCurSel())
	{
	case 8: //Apollo 14
		PIOSDataSetBox.SetWindowTextW(L"NBY71");
		break;
	}
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

	if (RopeNameBox.GetCurSel() <= CMC_COMANCHE072)
	{
		TRUNSFBox.SetReadOnly(true);
		SHAFTSFBox.SetReadOnly(true);
	}
	else
	{
		TRUNSFBox.SetReadOnly(false);
		SHAFTSFBox.SetReadOnly(false);
	}

	//TBD: Disable R2ModelBox somehow for Colossus 249 and earlier

	//SPS Performance

	//Colossus 237: None
	//Colossus 249: EMDOT
	//Comanche 45: EMDOT, EK1VAL, FANG
	//Comanche 55-108: EMDOT, EK1VAL, EK2VAL, EK3VAL, FANG
	//Artemis 72: EMDOT, EIMP1SEC, EFIMP01, EFIMP16

	switch (RopeNameBox.GetCurSel())
	{
	case CMC_COLOSSUS237:
	case CMC_COLOSSUS249:
		if (RopeNameBox.GetCurSel() == CMC_COLOSSUS237)
		{
			EMDOTBox.SetReadOnly(true);
		}
		else
		{
			EMDOTBox.SetReadOnly(false);
		}
		MinImp1Box.SetReadOnly(true);
		MinImp2Box.SetReadOnly(true);
		MinImp3Box.SetReadOnly(true);
		MinImp4Box.SetReadOnly(true);
		MinImp1Label.SetWindowTextW(L"NONE");
		MinImp2Label.SetWindowTextW(L"NONE");
		MinImp3Label.SetWindowTextW(L"NONE");
		MinImp4Label.SetWindowTextW(L"NONE");
		MinImp1Unit.SetWindowTextW(L"");
		MinImp2Unit.SetWindowTextW(L"");
		MinImp3Unit.SetWindowTextW(L"");
		MinImp4Unit.SetWindowTextW(L"");
		m_ToolTip.UpdateTipText(_T("None"), &MinImp1Box);
		m_ToolTip.UpdateTipText(_T("None"), &MinImp2Box);
		m_ToolTip.UpdateTipText(_T("None"), &MinImp3Box);
		m_ToolTip.UpdateTipText(_T("None"), &MinImp4Box);
		break;
	case CMC_COMANCHE045:
		EMDOTBox.SetReadOnly(false);
		MinImp1Box.SetReadOnly(false);
		MinImp2Box.SetReadOnly(true);
		MinImp3Box.SetReadOnly(true);
		MinImp4Box.SetReadOnly(false);
		MinImp1Label.SetWindowTextW(L"EK1VAL");
		MinImp2Label.SetWindowTextW(L"NONE");
		MinImp3Label.SetWindowTextW(L"NONE");
		MinImp4Label.SetWindowTextW(L"FANG");
		MinImp1Unit.SetWindowTextW(L"lbf/sec");
		MinImp2Unit.SetWindowTextW(L"");
		MinImp3Unit.SetWindowTextW(L"");
		MinImp4Unit.SetWindowTextW(L"lbf");
		m_ToolTip.UpdateTipText(_T("SPS impulse acquired from a one second burn"), &MinImp1Box);
		m_ToolTip.UpdateTipText(_T("None"), &MinImp2Box);
		m_ToolTip.UpdateTipText(_T("None"), &MinImp3Box);
		m_ToolTip.UpdateTipText(_T("SPS thrust used to estimate burn time when burn time is less than 6 sec"), &MinImp4Box);

		MinImp1Box.SetWindowTextW(L"19840.1");
		MinImp4Box.SetWindowTextW(L"20180.0");
		break;
	case CMC_COMANCHE055:
	case CMC_COMANCHE067:
	case CMC_COMANCHE072:
	case CMC_COMANCHE108:
		EMDOTBox.SetReadOnly(false);
		MinImp1Box.SetReadOnly(false);
		MinImp2Box.SetReadOnly(false);
		MinImp3Box.SetReadOnly(false);
		MinImp4Box.SetReadOnly(false);
		MinImp1Label.SetWindowTextW(L"EK1VAL");
		MinImp2Label.SetWindowTextW(L"EK2VAL");
		MinImp3Label.SetWindowTextW(L"EK3VAL");
		MinImp4Label.SetWindowTextW(L"FANG");
		MinImp1Unit.SetWindowTextW(L"lbf/sec");
		MinImp2Unit.SetWindowTextW(L"lbf/sec");
		MinImp3Unit.SetWindowTextW(L"lbf");
		MinImp4Unit.SetWindowTextW(L"lbf");
		m_ToolTip.UpdateTipText(_T("SPS impulse acquired from a one second burn"), &MinImp1Box);
		m_ToolTip.UpdateTipText(_T("SPS minimum impulse constant used to estimate burn time when burn time is less than 1.0 sec"), &MinImp2Box);
		m_ToolTip.UpdateTipText(_T("SPS minimum constant equal to the slope of the minimum impulse curve. Used to estimate burn time when burn time is less than 1.0 sec"), &MinImp3Box);
		m_ToolTip.UpdateTipText(_T("SPS thrust used to estimate burn time when burn time is less than 6 sec"), &MinImp4Box);

		MinImp1Box.SetWindowTextW(L"20143.2");
		MinImp2Box.SetWindowTextW(L"4909.1");
		MinImp3Box.SetWindowTextW(L"25454.5");
		MinImp4Box.SetWindowTextW(L"20240.0");
		break;
	default:
		EMDOTBox.SetReadOnly(false);
		MinImp1Box.SetReadOnly(false);
		MinImp2Box.SetReadOnly(false);
		MinImp3Box.SetReadOnly(false);
		MinImp4Box.SetReadOnly(true);
		MinImp1Label.SetWindowTextW(L"EIMP1SEC");
		MinImp2Label.SetWindowTextW(L"EFIMP01");
		MinImp3Label.SetWindowTextW(L"EFIMP16");
		MinImp4Label.SetWindowTextW(L"NONE");
		MinImp1Unit.SetWindowTextW(L"lbf/sec");
		MinImp2Unit.SetWindowTextW(L"lbf");
		MinImp3Unit.SetWindowTextW(L"lbf");
		MinImp4Unit.SetWindowTextW(L"");
		m_ToolTip.UpdateTipText(_T("Impulse from first second of SPS thrusting"), &MinImp1Box);
		m_ToolTip.UpdateTipText(_T("Slope of minimum impulse curve for SPS 0-1 second"), &MinImp2Box);
		m_ToolTip.UpdateTipText(_T("Slope of minimum impulse curve for SPS 1-6 seconds"), &MinImp3Box);
		m_ToolTip.UpdateTipText(_T("None"), &MinImp4Box);

		MinImp1Box.SetWindowTextW(L"19538.26");
		MinImp2Box.SetWindowTextW(L"24447.36");
		MinImp3Box.SetWindowTextW(L"19989.0");
		break;
	}

	//PIOS Data Set
	switch (RopeNameBox.GetCurSel())
	{
	case CMC_COLOSSUS237:
	case CMC_COLOSSUS249:
	case CMC_COMANCHE045:
		PIOSDataSetBox.SetWindowTextW(L"NBY1969");
		break;
	case CMC_COMANCHE055:
		PIOSDataSetBox.SetWindowTextW(L"NBY1970_V1");
		break;
	case CMC_COMANCHE067:
	case CMC_COMANCHE072:
		PIOSDataSetBox.SetWindowTextW(L"NBY1970_V2");
		break;
	case CMC_COMANCHE108:
		PIOSDataSetBox.SetWindowTextW(L"NBY1971");
		break;
	default:
		PIOSDataSetBox.SetWindowTextW(L"NBY1972");
		break;
	}
}

void CAGCPadloadGeneratorGUIDlg::Apollo7Padload()
{
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

	PBIASXBox.SetWindowText(L"0.24");		// cm/sec^2
	PIPASCFXBox.SetWindowText(L"-300.0");	// ppm
	PBIASYBox.SetWindowText(L"0.24");		// cm/sec^2
	PIPASCFYBox.SetWindowText(L"-190.0");	// ppm
	PBIASZBox.SetWindowText(L"0.14");		// cm/sec^2
	PIPASCFZBox.SetWindowText(L"-340.0");	// ppm
	NBDXBox.SetWindowText(L"-0.5");			// meru
	NBDYBox.SetWindowText(L"0.0");			// meru
	NBDZBox.SetWindowText(L"-0.6");			// meru
	ADIAXBox.SetWindowText(L"8.2");			// meru/g
	ADIAYBox.SetWindowText(L"11.6");		// meru/g
	ADIAZBox.SetWindowText(L"20.8");		// meru/g
	ADSRAXBox.SetWindowText(L"3.9");		// meru/g
	ADSRAYBox.SetWindowText(L"-0.4");		// meru/g
	ADSRAZBox.SetWindowText(L"-8.8");		// meru/g
}

void CAGCPadloadGeneratorGUIDlg::Apollo8Padload()
{
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

	PBIASXBox.SetWindowText(L"-0.0013");	// cm/sec^2
	PIPASCFXBox.SetWindowText(L"-76.57");	// ppm
	PBIASYBox.SetWindowText(L"0.803");		// cm/sec^2
	PIPASCFYBox.SetWindowText(L"-329.14");	// ppm
	PBIASZBox.SetWindowText(L"0.631");		// cm/sec^2
	PIPASCFZBox.SetWindowText(L"-200.71");	// ppm
	NBDXBox.SetWindowText(L"2.18");			// meru
	NBDYBox.SetWindowText(L"2.62");			// meru
	NBDZBox.SetWindowText(L"3.22");			// meru
	ADIAXBox.SetWindowText(L"17.88");		// meru/g
	ADIAYBox.SetWindowText(L"-3.22");		// meru/g
	ADIAZBox.SetWindowText(L"28.23");		// meru/g
	ADSRAXBox.SetWindowText(L"-2.24");		// meru/g
	ADSRAYBox.SetWindowText(L"1.66");		// meru/g
	ADSRAZBox.SetWindowText(L"-3.22");		// meru/g
}

void CAGCPadloadGeneratorGUIDlg::Apollo9Padload()
{
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
	EMDOTBox.SetWindowTextW(L"65.272");

	PBIASXBox.SetWindowText(L"0.64");		// cm/sec^2
	PIPASCFXBox.SetWindowText(L"-140.0");	// ppm
	PBIASYBox.SetWindowText(L"-0.1");		// cm/sec^2
	PIPASCFYBox.SetWindowText(L"-330.0");	// ppm
	PBIASZBox.SetWindowText(L"0.44");		// cm/sec^2
	PIPASCFZBox.SetWindowText(L"-280.0");	// ppm
	NBDXBox.SetWindowText(L"2.4");			// meru
	NBDYBox.SetWindowText(L"0.0");			// meru
	NBDZBox.SetWindowText(L"2.4");			// meru
	ADIAXBox.SetWindowText(L"5.0");			// meru/g
	ADIAYBox.SetWindowText(L"8.0");			// meru/g
	ADIAZBox.SetWindowText(L"-18.0");		// meru/g
	ADSRAXBox.SetWindowText(L"7.0");		// meru/g
	ADSRAYBox.SetWindowText(L"9.0");		// meru/g
	ADSRAZBox.SetWindowText(L"-4.0");		// meru/g
}

void CAGCPadloadGeneratorGUIDlg::Apollo10Padload()
{
	EphemerisSpanBox.SetWindowTextW(L"10.5");
	LaunchMJDInput.SetWindowTextW(L"40359.7006944");
	Launchpad.SetCurSel(2); //LC-39B
	RTEDBox.SetWindowTextW(L"1.6602637");
	LSAltitudeBox.SetWindowTextW(L"-3073.26");
	LSLatitudeBox.SetWindowTextW(L"0.71388");
	LSLongitudeBox.SetWindowTextW(L"23.707773");
	EMSAltBox.SetWindowTextW(L"300997.3");
	LaunchAzimuthBox.SetWindowTextW(L"72.028");
	HORIZALTBox.SetWindowTextW(L"24000");
	ALTVARBox.SetWindowText(L"1.5258e-05");
	CSMMASSBox.SetWindowTextW(L"63904.0");
	LEMMASSBox.SetWindowTextW(L"31579.7");
	PACTOFFBox.SetWindowTextW(L"-1.48");
	YACTOFFBox.SetWindowTextW(L"1.35");
	EMDOTBox.SetWindowTextW(L"66.957");
	MinImp1Box.SetWindowTextW(L"19840.1");
	MinImp4Box.SetWindowTextW(L"20180.0");
	LADPADBox.SetWindowTextW(L"0.27");
	LODPADBox.SetWindowTextW(L"0.207");
	ALFAPADBox.SetWindowTextW(L"-19.6");
	P37RANGEBox.SetWindowTextW(L"1236.0");
	RVARMINBox.SetWindowText(L"900");
	WMIDPOSBox.SetWindowText(L"5700");
	WMIDVELBox.SetWindowText(L"5.73"); //TBD: Gives the right octal, but improve this

	PBIASXBox.SetWindowText(L"-0.27");		// cm/sec^2
	PIPASCFXBox.SetWindowText(L"-100.0");	// ppm
	PBIASYBox.SetWindowText(L"-0.07");		// cm/sec^2
	PIPASCFYBox.SetWindowText(L"-230.0");	// ppm
	PBIASZBox.SetWindowText(L"-0.05");		// cm/sec^2
	PIPASCFZBox.SetWindowText(L"-80.0");	// ppm
	NBDXBox.SetWindowText(L"0.4");			// meru
	NBDYBox.SetWindowText(L"-1.3");			// meru
	NBDZBox.SetWindowText(L"1.2");			// meru
	ADIAXBox.SetWindowText(L"1.0");			// meru/g
	ADIAYBox.SetWindowText(L"13.0");		// meru/g
	ADIAZBox.SetWindowText(L"11.0");		// meru/g
	ADSRAXBox.SetWindowText(L"10.0");		// meru/g
	ADSRAYBox.SetWindowText(L"3.0");		// meru/g
	ADSRAZBox.SetWindowText(L"7.0");		// meru/g
}

void CAGCPadloadGeneratorGUIDlg::Apollo11Padload(int LaunchDay)
{
	if (LaunchDay == 16)
	{
		LaunchMJDInput.SetWindowTextW(L"40418.5638889");
		LSAltitudeBox.SetWindowTextW(L"-2671.9684");
		LSLatitudeBox.SetWindowTextW(L"0.691395");
		LSLongitudeBox.SetWindowTextW(L"23.71689");
		LaunchAzimuthBox.SetWindowTextW(L"72.05897");
	}
	else if (LaunchDay == 18)
	{
		LaunchMJDInput.SetWindowTextW(L"40420.647222");
		LSAltitudeBox.SetWindowTextW(L"-1817.3");
		LSLatitudeBox.SetWindowTextW(L"0.35277778");
		LSLongitudeBox.SetWindowTextW(L"-1.29916667");
		LaunchAzimuthBox.SetWindowTextW(L"89.295");
	}
	else
	{
		LaunchMJDInput.SetWindowTextW(L"40423.672905");
		LSAltitudeBox.SetWindowTextW(L"-2314.0");
		LSLatitudeBox.SetWindowTextW(L"1.67805556");
		LSLongitudeBox.SetWindowTextW(L"-41.89916667");
		LaunchAzimuthBox.SetWindowTextW(L"94.6775");
	}

	EphemerisSpanBox.SetWindowTextW(L"10.5");

	Launchpad.SetCurSel(1); //LC-39A
	RTEDBox.SetWindowTextW(L"1.6602637");
	EMSAltBox.SetWindowTextW(L"294084.3");
	HORIZALTBox.SetWindowTextW(L"24000");
	ALTVARBox.SetWindowTextW(L"1.5258e-05");
	CSMMASSBox.SetWindowTextW(L"63386.7");
	LEMMASSBox.SetWindowTextW(L"33275.6");
	PACTOFFBox.SetWindowTextW(L"-1.541");
	YACTOFFBox.SetWindowTextW(L"1.321");
	LADPADBox.SetWindowTextW(L"0.27");
	LODPADBox.SetWindowTextW(L"0.207");
	ALFAPADBox.SetWindowTextW(L"-19.55");
	P37RANGEBox.SetWindowTextW(L"1221.5");
	EMDOTBox.SetWindowTextW(L"66.64");
	MinImp1Box.SetWindowTextW(L"20143.2");
	MinImp2Box.SetWindowTextW(L"4909.1");
	MinImp3Box.SetWindowTextW(L"25454.5");
	MinImp4Box.SetWindowTextW(L"20240.0");

	PBIASXBox.SetWindowText(L"-0.26");		// cm/sec^2
	PIPASCFXBox.SetWindowText(L"40.0");		// ppm
	PBIASYBox.SetWindowText(L"-0.13");		// cm/sec^2
	PIPASCFYBox.SetWindowText(L"-80.0");	// ppm
	PBIASZBox.SetWindowText(L"0.14");		// cm/sec^2
	PIPASCFZBox.SetWindowText(L"-30.0");	// ppm
	NBDXBox.SetWindowText(L"-1.8");			// meru
	NBDYBox.SetWindowText(L"-0.6");			// meru
	NBDZBox.SetWindowText(L"-0.2");			// meru
	ADIAXBox.SetWindowText(L"15.0");		// meru/g
	ADIAYBox.SetWindowText(L"5.0");			// meru/g
	ADIAZBox.SetWindowText(L"1.0");			// meru/g
	ADSRAXBox.SetWindowText(L"-6.0");		// meru/g
	ADSRAYBox.SetWindowText(L"3.0");		// meru/g
	ADSRAZBox.SetWindowText(L"5.0");		// meru/g

	agc.BLOCKII.POLYNUM[0] = 2.667379e-1;
	agc.BLOCKII.POLYNUM[1] = 0.3714781e-1;
	agc.BLOCKII.POLYNUM[2] = 0.917483e-2;
	agc.BLOCKII.POLYNUM[3] = 0.7781624e-4;
	agc.BLOCKII.POLYNUM[4] = -0.2593124e-5;
	agc.BLOCKII.POLYNUM[5] = 0.1818239e-7;
	agc.BLOCKII.POLYNUM[6] = -0.4160838e-10;
	agc.BLOCKII.RPSTART = 11.85;
	agc.BLOCKII.POLYSTOP = 147.25;
}

void CAGCPadloadGeneratorGUIDlg::Apollo12Padload()
{
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
	EMDOTBox.SetWindowTextW(L"66.95");
	MinImp1Box.SetWindowTextW(L"19293.0");
	MinImp2Box.SetWindowTextW(L"4909.1");
	MinImp3Box.SetWindowTextW(L"25454.5");
	MinImp4Box.SetWindowTextW(L"20260.0");

	PBIASXBox.SetWindowText(L"-0.09");		// cm/sec^2
	PIPASCFXBox.SetWindowText(L"-220.0");	// ppm
	PBIASYBox.SetWindowText(L"-0.09");		// cm/sec^2
	PIPASCFYBox.SetWindowText(L"-350.0");	// ppm
	PBIASZBox.SetWindowText(L"-0.16");		// cm/sec^2
	PIPASCFZBox.SetWindowText(L"-370.0");	// ppm
	NBDXBox.SetWindowText(L"-0.1");			// meru
	NBDYBox.SetWindowText(L"-0.1");			// meru
	NBDZBox.SetWindowText(L"0.1");			// meru
	ADIAXBox.SetWindowText(L"13.0");		// meru/g
	ADIAYBox.SetWindowText(L"0.0");			// meru/g
	ADIAZBox.SetWindowText(L"-1.0");		// meru/g
	ADSRAXBox.SetWindowText(L"-4.0");		// meru/g
	ADSRAYBox.SetWindowText(L"-4.0");		// meru/g
	ADSRAZBox.SetWindowText(L"-6.0");		// meru/g

	agc.BLOCKII.POLYNUM[0] = 1.1940272e-1;
	agc.BLOCKII.POLYNUM[1] = 8.0761057e-2;
	agc.BLOCKII.POLYNUM[2] = 1.4609205e-2;
	agc.BLOCKII.POLYNUM[3] = -1.9841095e-4;
	agc.BLOCKII.POLYNUM[4] = 1.7966877e-6;
	agc.BLOCKII.POLYNUM[5] = -1.1164227e-8;
	agc.BLOCKII.POLYNUM[6] = 2.9494747e-11;
	agc.BLOCKII.RPSTART = 12.6;
	agc.BLOCKII.POLYSTOP = 145.4;
}

void CAGCPadloadGeneratorGUIDlg::Apollo13Padload()
{
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
	EMDOTBox.SetWindowTextW(L"67.6");
	MinImp1Box.SetWindowTextW(L"18848.0");
	MinImp2Box.SetWindowTextW(L"6606.5");
	MinImp3Box.SetWindowTextW(L"25454.5");
	MinImp4Box.SetWindowTextW(L"20385.0");

	agc.BLOCKII.POLYNUM[0] = -7.646894e-2;
	agc.BLOCKII.POLYNUM[1] = 1.79988e-1;
	agc.BLOCKII.POLYNUM[2] = 9.298907e-3;
	agc.BLOCKII.POLYNUM[3] = -1.49242e-4;
	agc.BLOCKII.POLYNUM[4] = 2.078269e-6;
	agc.BLOCKII.POLYNUM[5] = -1.601873e-8;
	agc.BLOCKII.POLYNUM[6] = 4.401981e-11;
	agc.BLOCKII.RPSTART = 11.85;
	agc.BLOCKII.POLYSTOP = 149.5;
}

void CAGCPadloadGeneratorGUIDlg::Apollo14Padload()
{
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
	EMDOTBox.SetWindowTextW(L"67.5");
	MinImp1Box.SetWindowTextW(L"19965.7");
	MinImp2Box.SetWindowTextW(L"4909.1");
	MinImp3Box.SetWindowTextW(L"24874.8");
	MinImp4Box.SetWindowTextW(L"20390.0");
	TRUNSFBox.SetWindowTextW(L"1269760.0");
	SHAFTSFBox.SetWindowTextW(L"659456.0");

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
	EphemerisSpanBox.SetWindowTextW(L"14.5");
	LaunchMJDInput.SetWindowTextW(L"41158.565277778");
	Launchpad.SetCurSel(1); //LC-39A
	RTEDBox.SetWindowTextW(L"1.6602637");
	LSAltitudeBox.SetWindowTextW(L"-3550.0");
	LSLatitudeBox.SetWindowTextW(L"26.07389");
	LSLongitudeBox.SetWindowTextW(L"3.65389");
	EMSAltBox.SetWindowTextW(L"300781.4");
	LaunchAzimuthBox.SetWindowTextW(L"80.08868");
	LATSPLBox.SetWindowTextW(L"20.3");
	LNGSPLBox.SetWindowTextW(L"-19.5");
	HORIZALTBox.SetWindowTextW(L"28000");
	CSMMASSBox.SetWindowTextW(L"66915.2");
	LEMMASSBox.SetWindowTextW(L"36206.2");
	PACTOFFBox.SetWindowTextW(L"-0.524");
	YACTOFFBox.SetWindowTextW(L"1.895");
	LADPADBox.SetWindowTextW(L"0.27");
	LODPADBox.SetWindowTextW(L"0.207");
	ALFAPADBox.SetWindowTextW(L"-18.64");
	P37RANGEBox.SetWindowTextW(L"1114.55");
	EMDOTBox.SetWindowTextW(L"67.5");
	MinImp1Box.SetWindowTextW(L"19538.26");
	MinImp2Box.SetWindowTextW(L"24447.36");
	MinImp3Box.SetWindowTextW(L"19989.0");
	TRUNSFBox.SetWindowTextW(L"1269760.0");
	SHAFTSFBox.SetWindowTextW(L"647168.0");

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
	LATSPLBox.SetWindowTextW(L"26.5");
	LNGSPLBox.SetWindowTextW(L"-17.0");
	EMDOTBox.SetWindowTextW(L"65.97842");
	MinImp1Box.SetWindowTextW(L"19286.838");
	MinImp2Box.SetWindowTextW(L"24195.938");
	MinImp3Box.SetWindowTextW(L"20779.0");
	TRUNSFBox.SetWindowTextW(L"1269760.0");
	SHAFTSFBox.SetWindowTextW(L"651264.0");

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
	LATSPLBox.SetWindowTextW(L"26.35");
	LNGSPLBox.SetWindowTextW(L"-17.1");
	EMDOTBox.SetWindowTextW(L"66.23333");
	MinImp1Box.SetWindowTextW(L"19575.168");
	MinImp2Box.SetWindowTextW(L"24484.268");
	MinImp3Box.SetWindowTextW(L"20281.0");
	TRUNSFBox.SetWindowTextW(L"1286144.0");
	SHAFTSFBox.SetWindowTextW(L"659456.0");

	agc.BLOCKII.POLYNUM[0] = -3.7352617e-1;
	agc.BLOCKII.POLYNUM[1] = 2.4436145e-1;
	agc.BLOCKII.POLYNUM[2] = 6.2373017e-3;
	agc.BLOCKII.POLYNUM[3] = -9.0587464e-5;
	agc.BLOCKII.POLYNUM[4] = 1.685017e-6;
	agc.BLOCKII.POLYNUM[5] = -1.6095034e-8;
	agc.BLOCKII.POLYNUM[6] = 4.9817897e-11;
	agc.BLOCKII.RPSTART = 12.25;
	agc.BLOCKII.POLYSTOP = 147.75;

	/*agc.CMCDATA.IMUBiasCompensation.PBIASX = -0.19; // cm/sec^2
	agc.CMCDATA.IMUBiasCompensation.PIPASCFX = -750.0; // ppm
	agc.CMCDATA.IMUBiasCompensation.PBIASY = -0.07; // cm/sec^2
	agc.CMCDATA.IMUBiasCompensation.PIPASCFY = -210.0; // ppm
	agc.CMCDATA.IMUBiasCompensation.PBIASZ = 0.71; // cm/sec^2
	agc.CMCDATA.IMUBiasCompensation.PIPASCFZ = -810.0; // ppm
	agc.CMCDATA.IMUBiasCompensation.NBDX = 0.6; // meru
	agc.CMCDATA.IMUBiasCompensation.NBDY = -0.1; // meru
	agc.CMCDATA.IMUBiasCompensation.NBDZ = 1.2; // meru
	agc.CMCDATA.IMUBiasCompensation.ADIAX = 7.0; // meru/g
	agc.CMCDATA.IMUBiasCompensation.ADIAY = 13.0; // meru/g
	agc.CMCDATA.IMUBiasCompensation.ADIAZ = -6.0; // meru/g
	agc.CMCDATA.IMUBiasCompensation.ADSRAX = -6.0; // meru/g
	agc.CMCDATA.IMUBiasCompensation.ADSRAY = 6.0; // meru/g
	agc.CMCDATA.IMUBiasCompensation.ADSRAZ = -4.0; // meru/g*/
}