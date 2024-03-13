#pragma once

#include "AGCPadloadGenerator.h"

// LGCPadloadGenerator-Dialog

class LGCPadloadGenerator : public CDialogEx
{
	DECLARE_DYNAMIC(LGCPadloadGenerator)

public:
	LGCPadloadGenerator(CWnd* pParent = nullptr);   // Standardkonstruktor
	virtual ~LGCPadloadGenerator();

	BOOL PreTranslateMessage(MSG* pMsg);

// Dialogfelddaten
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_DIALOG2 };
#endif

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV-Unterstützung

	virtual BOOL OnInitDialog();
	DECLARE_MESSAGE_MAP()

	void UpdateTotalMass();
	void UpdateRopeSpecificEditFields();

	CToolTipCtrl m_ToolTip;
public:

	AGCPadloadGenerator agc;

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
	CEdit LMMassBox;
	CEdit CSMMassBox;
	CEdit TotalMassBox;
	CButton DockedBox;
	afx_msg void OnEnChangeEdit6();
	afx_msg void OnEnChangeEdit7();
	afx_msg void OnBnClickedCheck1();
	CEdit WRENDPOSBox;
	CEdit WRENDVELBox;
	CEdit WSHAFTBox;
	CEdit WTRUNBox;
	CEdit RMAXBox;
	CEdit VMAXBox;
	CEdit SHAFTVARBox;
	CEdit TRUNVARBox;
	CEdit WSURFPOSBox;
	CEdit WSURFVELBox;
	CEdit HIASCENTBox;
	CEdit AGSKBox;
	CEdit ROLLTIMEBox;
	CEdit PITCHTIMEBox;
	CEdit TLANDBox;
	CEdit ABSC0Box;
	CEdit ABSC1Box;
	CEdit ABSC2Box;
	CEdit ABSC3Box;
	CEdit ABSC4Box;
	CEdit SLOPE0Box;
	CEdit SLOPE1Box;
	CEdit SLOPE2Box;
	CEdit SLOPE3Box;
	CEdit SLOPE4Box;
	CEdit IGNAOSQBox;
	CEdit IGNAOSRBox;
	CEdit VIGNBox;
	CEdit RIGNXBox;
	CEdit RIGNZBox;
	CEdit KIGNXBox;
	CEdit KIGNYBox;
	CEdit KIGNVBox;
	CEdit J1PARMBox;
	CEdit K1PARMBox;
	CEdit J2PARMBox;
	CEdit K2PARMBox;
	CEdit THETCRITBox;
	CEdit RAMINBox;
	CEdit DELTTFAPBox;
	CEdit OutputBox;
	CEdit PIOSDataSetBox;
};
