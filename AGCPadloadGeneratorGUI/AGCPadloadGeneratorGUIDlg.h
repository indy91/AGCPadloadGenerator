
// AGCPadloadGeneratorGUIDlg.h: Headerdatei
//

#pragma once

#include "AGCPadloadGenerator.h"

// CAGCPadloadGeneratorGUIDlg-Dialogfeld
class CAGCPadloadGeneratorGUIDlg : public CDialogEx
{
// Konstruktion
public:
	CAGCPadloadGeneratorGUIDlg(CWnd* pParent = nullptr);	// Standardkonstruktor

// Dialogfelddaten
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_AGCPADLOADGENERATORGUI_DIALOG };
#endif

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV-Unterstützung


// Implementierung
protected:
	HICON m_hIcon;

	// Generierte Funktionen für die Meldungstabellen
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()

	void String2Text(std::string val, CEdit *ed);
	void String2Text(std::string val, CComboBox *ed);
public:
	CEdit LaunchMJDInput;
	afx_msg void OnBnClickedOk();
	afx_msg void OnCbnSelchangeCombo1();
	afx_msg void OnCbnSelchangeCombo2();
	afx_msg void OnCbnSelchangeCombo3();

	AGCPadloadGenerator agc;

	CComboBox Launchpad;
	CString LaunchPadValue;
	CEdit TLANDBox;
	CEdit LSLatitudeBox;
	CEdit LSLongitudeBox;
	CEdit LSAltitudeBox;
	CEdit RTEDBox;
	CComboBox RopeNameBox;
	CEdit EMSAltBox;
	CEdit UNITWBox;
	CEdit EphemerisSpanBox;
	CComboBox MissionBox;
	CString MissionNameValue;
	CEdit LaunchAzimuthBox;
	CEdit CDUCHKWDBox;
	CEdit HORIZALTBox;
	CEdit ALTVARBox;
	CEdit WRENDPOSBox;
	CEdit WRENDVELBox;
	CEdit RMAXBox;
	CEdit VMAXBox;
};
