#pragma once

#include "AGCPadloadGenerator.h"

// BlockIPadloadGenerator dialog

class BlockIPadloadGenerator : public CDialogEx
{
	DECLARE_DYNAMIC(BlockIPadloadGenerator)

public:
	BlockIPadloadGenerator(CWnd* pParent = nullptr);   // standard constructor
	virtual ~BlockIPadloadGenerator();

// Dialog Data
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_DIALOG3 };
#endif

protected:
	virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	AGCPadloadGenerator agc;

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedOk();
	CComboBox Launchpad;
	CComboBox MissionBox;
	CComboBox RopeNameBox;
	afx_msg void OnCbnSelchangeCombo3();
	CEdit LaunchMJDInput;
	CEdit LaunchAzimuthBox;
	CEdit SPS1EccentricityBox;
	CEdit SPS2EccentricityBox;
	CEdit SPS1SMABox;
	CEdit SPS2SMABox;
	CEdit TAtlanticBox;
	CEdit TPacificBox;
	CEdit AtlanticLatitudeBox;
	CEdit AtlanticLongitudeBox;
	CEdit PacificLatitudeBox;
	CEdit PacificLongitudeBox;
};
