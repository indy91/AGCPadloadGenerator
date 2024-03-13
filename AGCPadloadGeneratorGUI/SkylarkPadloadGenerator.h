#pragma once

#include "AGCPadloadGenerator.h"

// SkylarkPadloadGenerator dialog

class SkylarkPadloadGenerator : public CDialogEx
{
	DECLARE_DYNAMIC(SkylarkPadloadGenerator)

public:
	SkylarkPadloadGenerator(CWnd* pParent = nullptr);   // standard constructor
	virtual ~SkylarkPadloadGenerator();

// Dialog Data
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_DIALOG4 };
#endif

protected:
	virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()

	AGCPadloadGenerator agc;
public:
	afx_msg void OnBnClickedOk();
	CEdit LaunchMJDInput;
	CComboBox MissionBox;
	afx_msg void OnCbnSelchangeCombo1();
	CEdit LaunchAzimuthBox;
	CEdit CSMMASSBox;
	CEdit LEMMASSBox;
	CEdit CDUCHKWDBox;
	CEdit EMSALTBox;
	CEdit LADPADBox;
	CEdit LODPADBox;
	CEdit ALFAPADBox;
	CEdit TTPIBox;
	CEdit PACTOFFBox;
	CEdit YACTOFFBox;
	CEdit LATSPLBox;
	CEdit LNGSPLBox;
	CEdit C12ALPHABox;
	CEdit ECPBox;
	CEdit ECYWBox;
	CEdit ALPHAPBox;
	CEdit ALPHAYWBox;
	CEdit KMJDCKDBox;
	CEdit KMJ1DCKDBox;
	CEdit KMJ2DCKDBox;
	CEdit JMDCKDBox;
	CEdit JM1DCKDBox;
	CEdit JM2DCKDBox;
	CEdit CH6FAILBox;
	CEdit DKRATEBox;
	CEdit OutputBox;
};
