// SkylarkPadloadGenerator.cpp : implementation file
//

#include "pch.h"
#include "AGCPadloadGeneratorGUI.h"
#include "SkylarkPadloadGenerator.h"
#include "afxdialogex.h"


// SkylarkPadloadGenerator dialog

IMPLEMENT_DYNAMIC(SkylarkPadloadGenerator, CDialogEx)

SkylarkPadloadGenerator::SkylarkPadloadGenerator(CWnd* pParent /*=nullptr*/)
	: CDialogEx(IDD_DIALOG4, pParent)
{

}

SkylarkPadloadGenerator::~SkylarkPadloadGenerator()
{
}

void SkylarkPadloadGenerator::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_EDIT1, LaunchMJDInput);
	DDX_Control(pDX, IDC_COMBO1, MissionBox);
	DDX_Control(pDX, IDC_EDIT2, LaunchAzimuthBox);
	DDX_Control(pDX, IDC_EDIT3, CSMMASSBox);
	DDX_Control(pDX, IDC_EDIT4, LEMMASSBox);
	DDX_Control(pDX, IDC_EDIT11, CDUCHKWDBox);
	DDX_Control(pDX, IDC_EDIT7, EMSALTBox);
	DDX_Control(pDX, IDC_EDIT23, LADPADBox);
	DDX_Control(pDX, IDC_EDIT24, LODPADBox);
	DDX_Control(pDX, IDC_EDIT25, ALFAPADBox);
	DDX_Control(pDX, IDC_EDIT5, TTPIBox);
	DDX_Control(pDX, IDC_EDIT22, PACTOFFBox);
	DDX_Control(pDX, IDC_EDIT52, YACTOFFBox);
	DDX_Control(pDX, IDC_EDIT18, LATSPLBox);
	DDX_Control(pDX, IDC_EDIT19, LNGSPLBox);
	DDX_Control(pDX, IDC_EDIT6, C12ALPHABox);
	DDX_Control(pDX, IDC_EDIT53, ECPBox);
	DDX_Control(pDX, IDC_EDIT54, ECYWBox);
	DDX_Control(pDX, IDC_EDIT55, ALPHAPBox);
	DDX_Control(pDX, IDC_EDIT56, ALPHAYWBox);
	DDX_Control(pDX, IDC_EDIT57, KMJDCKDBox);
	DDX_Control(pDX, IDC_EDIT58, KMJ1DCKDBox);
	DDX_Control(pDX, IDC_EDIT59, KMJ2DCKDBox);
	DDX_Control(pDX, IDC_EDIT60, JMDCKDBox);
	DDX_Control(pDX, IDC_EDIT61, JM1DCKDBox);
	DDX_Control(pDX, IDC_EDIT62, JM2DCKDBox);
	DDX_Control(pDX, IDC_EDIT63, CH6FAILBox);
	DDX_Control(pDX, IDC_EDIT64, DKRATEBox);
}


BEGIN_MESSAGE_MAP(SkylarkPadloadGenerator, CDialogEx)
	ON_BN_CLICKED(IDOK, &SkylarkPadloadGenerator::OnBnClickedOk)
	ON_CBN_SELCHANGE(IDC_COMBO1, &SkylarkPadloadGenerator::OnCbnSelchangeCombo1)
END_MESSAGE_MAP()

BOOL SkylarkPadloadGenerator::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	LaunchMJDInput.SetWindowText(L"0.0");

	MissionBox.AddString(L"Manual");
	MissionBox.AddString(L"Skylab 2");
	MissionBox.AddString(L"Skylab 3");
	MissionBox.AddString(L"Skylab 4");
	MissionBox.AddString(L"ASTP");
	MissionBox.SetCurSel(0);

	LaunchAzimuthBox.SetWindowText(L"72.0");

	CSMMASSBox.SetWindowText(L"30000.0");
	LEMMASSBox.SetWindowText(L"0.0");

	CDUCHKWDBox.SetWindowText(L"5");

	LADPADBox.SetWindowText(L"0.3");
	LODPADBox.SetWindowText(L"0.18");
	ALFAPADBox.SetWindowText(L"-20.0");
	PACTOFFBox.SetWindowText(L"0.0");
	YACTOFFBox.SetWindowText(L"0.0");
	LATSPLBox.SetWindowText(L"0.0");
	LNGSPLBox.SetWindowText(L"0.0");

	EMSALTBox.SetWindowText(L"290000.0");
	TTPIBox.SetWindowText(L"0.0");

	C12ALPHABox.SetWindowText(L"6.3454574");
	ECPBox.SetWindowText(L"0.164185153");
	ECYWBox.SetWindowText(L"-0.162029721");
	ALPHAPBox.SetWindowText(L"-0.283407807");
	ALPHAYWBox.SetWindowText(L"0.647582933");
	KMJDCKDBox.SetWindowText(L"0.054069332");
	KMJ1DCKDBox.SetWindowText(L"0.053174108");
	KMJ2DCKDBox.SetWindowText(L"0.055783856");
	JMDCKDBox.SetWindowText(L"18.86467");
	JM1DCKDBox.SetWindowText(L"116.83315");
	JM2DCKDBox.SetWindowText(L"112.8488");
	CH6FAILBox.SetWindowText(L"030");
	DKRATEBox.SetWindowText(L"0.2");

	return TRUE;
}
// SkylarkPadloadGenerator message handlers


void SkylarkPadloadGenerator::OnBnClickedOk()
{
	agc.RopeName = "Skylark048";
	agc.Pad = "LC-39B";

	//TBD
	agc.BLOCKII.POLYNUM[0] = 5.463972e-2;
	agc.BLOCKII.POLYNUM[1] = 2.216763e-1;
	agc.BLOCKII.POLYNUM[2] = 6.4867174e-3;
	agc.BLOCKII.POLYNUM[3] = 8.9707548e-6;
	agc.BLOCKII.POLYNUM[4] = -9.5787142e-7;
	agc.BLOCKII.POLYNUM[5] = 6.60410301e-9;
	agc.BLOCKII.POLYNUM[6] = -1.2743408e-11;
	agc.BLOCKII.RPSTART = 10.3;
	agc.BLOCKII.POLYSTOP = 120.75;

	agc.LaunchMJD = Utilities::Text2Double(&LaunchMJDInput);
	agc.LaunchAzimuth = Utilities::Text2Double(&LaunchAzimuthBox);
	agc.CDUCHKWD = Utilities::Text2Double(&CDUCHKWDBox) / 100.0;
	agc.BLOCKII.LAT_SPL = Utilities::Text2Double(&LATSPLBox);
	agc.BLOCKII.LNG_SPL = Utilities::Text2Double(&LNGSPLBox);
	agc.BLOCKII.CSMMass = Utilities::Text2Double(&CSMMASSBox);
	agc.BLOCKII.LMMass = Utilities::Text2Double(&LEMMASSBox);
	agc.BLOCKII.PACTOFF = Utilities::Text2Double(&PACTOFFBox);
	agc.BLOCKII.YACTOFF = Utilities::Text2Double(&YACTOFFBox);
	agc.BLOCKII.LADPAD = Utilities::Text2Double(&LADPADBox);
	agc.BLOCKII.LODPAD = Utilities::Text2Double(&LODPADBox);
	agc.BLOCKII.ALFAPAD = Utilities::Text2Double(&ALFAPADBox);
	agc.CMCDATA.TTPI = Utilities::Text2Double(&TTPIBox);
	agc.EMSALT = Utilities::Text2Double(&EMSALTBox);
	agc.CMCDATA.C12ALPHA = Utilities::Text2Double(&C12ALPHABox);
	agc.CMCDATA.ECP = Utilities::Text2Double(&ECPBox);
	agc.CMCDATA.ECYW = Utilities::Text2Double(&ECYWBox);
	agc.CMCDATA.ALPHAP = Utilities::Text2Double(&ALPHAPBox);
	agc.CMCDATA.ALPHAYW = Utilities::Text2Double(&ALPHAYWBox);
	agc.CMCDATA.KMJDCKD = Utilities::Text2Double(&KMJDCKDBox);
	agc.CMCDATA.KMJ1DCKD = Utilities::Text2Double(&KMJ1DCKDBox);
	agc.CMCDATA.KMJ2DCKD = Utilities::Text2Double(&KMJ2DCKDBox);
	agc.CMCDATA.JMDCKD = Utilities::Text2Double(&JMDCKDBox);
	agc.CMCDATA.JM1DCKD = Utilities::Text2Double(&JM1DCKDBox);
	agc.CMCDATA.JM2DCKD = Utilities::Text2Double(&JM2DCKDBox);
	agc.CMCDATA.CH6FAIL = Utilities::Text2Octal(&CH6FAILBox);
	agc.CMCDATA.DKRATE = Utilities::Text2Double(&DKRATEBox);

	agc.RunCMC();
}


void SkylarkPadloadGenerator::OnCbnSelchangeCombo1()
{
	switch (MissionBox.GetCurSel())
	{
	case 1: //Skylab 2
		//LaunchMJDInput.SetWindowTextW(L"41817.70805555555556"); //Planned
		LaunchMJDInput.SetWindowTextW(L"41827.54166667"); //Actual
		LaunchAzimuthBox.SetWindowText(L"47.035");
		CSMMASSBox.SetWindowText(L"30442.5");
		LEMMASSBox.SetWindowText(L"0.0");
		LADPADBox.SetWindowText(L"0.25");
		LODPADBox.SetWindowText(L"0.225");
		ALFAPADBox.SetWindowText(L"-20.49");
		LATSPLBox.SetWindowText(L"49.2");
		LNGSPLBox.SetWindowText(L"-7.3");
		TTPIBox.SetWindowText(L"24578.7");
		PACTOFFBox.SetWindowText(L"1.088");
		YACTOFFBox.SetWindowText(L"0.441");

		C12ALPHABox.SetWindowText(L"6.3454574");
		ECPBox.SetWindowText(L"0.164185153");
		ECYWBox.SetWindowText(L"-0.162029721");
		ALPHAPBox.SetWindowText(L"-0.283407807");
		ALPHAYWBox.SetWindowText(L"0.647582933");
		KMJDCKDBox.SetWindowText(L"0.054069332");
		KMJ1DCKDBox.SetWindowText(L"0.053174108");
		KMJ2DCKDBox.SetWindowText(L"0.055783856");
		JMDCKDBox.SetWindowText(L"18.86467");
		JM1DCKDBox.SetWindowText(L"116.83315");
		JM2DCKDBox.SetWindowText(L"112.8488");
		CH6FAILBox.SetWindowText(L"030");
		DKRATEBox.SetWindowText(L"0.2");
		break;
	case 2: //Skylab 3
		break;
	case 3: //Skylab 4
		break;
	case 4: //ASTP
		LaunchMJDInput.SetWindowTextW(L"42608.826389");
		LaunchAzimuthBox.SetWindowText(L"45.159");
		CSMMASSBox.SetWindowText(L"27885.0");
		LEMMASSBox.SetWindowText(L"0.0");
		LADPADBox.SetWindowText(L"0.25");
		LODPADBox.SetWindowText(L"0.225");
		ALFAPADBox.SetWindowText(L"-18.77");
		LATSPLBox.SetWindowText(L"50.2");
		LNGSPLBox.SetWindowText(L"-9.7");
		TTPIBox.SetWindowText(L"156265.0");
		PACTOFFBox.SetWindowText(L"0.7008");
		YACTOFFBox.SetWindowText(L"-0.6045");

		C12ALPHABox.SetWindowText(L"0.4278333333");
		ECPBox.SetWindowText(L"0.484558");
		ECYWBox.SetWindowText(L"-0.483093");
		ALPHAPBox.SetWindowText(L"0.0");
		ALPHAYWBox.SetWindowText(L"0.0");
		KMJDCKDBox.SetWindowText(L"2.000376");
		KMJ1DCKDBox.SetWindowText(L"0.19044");
		KMJ2DCKDBox.SetWindowText(L"0.19044");
		JMDCKDBox.SetWindowText(L"0.5006222222");
		JM1DCKDBox.SetWindowText(L"11.0592");
		JM2DCKDBox.SetWindowText(L"11.0592");
		CH6FAILBox.SetWindowText(L"030");
		DKRATEBox.SetWindowText(L"0.5");
		break;
	}
}
