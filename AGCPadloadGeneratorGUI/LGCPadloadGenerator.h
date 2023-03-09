#pragma once

#include "AGCPadloadGenerator.h"

// LGCPadloadGenerator-Dialog

class LGCPadloadGenerator : public CDialogEx
{
	DECLARE_DYNAMIC(LGCPadloadGenerator)

public:
	LGCPadloadGenerator(CWnd* pParent = nullptr);   // Standardkonstruktor
	virtual ~LGCPadloadGenerator();

// Dialogfelddaten
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_DIALOG2 };
#endif

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV-Unterstützung

	virtual BOOL OnInitDialog();
	DECLARE_MESSAGE_MAP()
public:

	AGCPadloadGenerator agc;

	CString RopeNameValue;
	CString MissionNameValue;

	CComboBox RopeNameBox;
	CComboBox MissionBox;

	afx_msg void OnBnClickedOk();
	afx_msg void OnCbnSelchangeCombo1();
	afx_msg void OnCbnSelchangeCombo2();
	CEdit LaunchMJDInput;
	CEdit EphemerisSpanBox;
	CEdit T504LMBox;
	CEdit UNITWBox;
	CEdit LSLatitudeBox;
	CEdit LSLongitudeBox;
	CEdit LSAltitudeBox;
	afx_msg void OnBnClickedCancel();
};
