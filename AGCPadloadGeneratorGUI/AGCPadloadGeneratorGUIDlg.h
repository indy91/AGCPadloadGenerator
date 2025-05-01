
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

	BOOL PreTranslateMessage(MSG* pMsg);

// Implementierung
protected:
	HICON m_hIcon;
	CToolTipCtrl m_ToolTip;

	// Generierte Funktionen für die Meldungstabellen
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()

	void String2Text(std::string val, CEdit *ed);
	void String2Text(std::string val, CComboBox *ed);

	void UpdateRopeSpecificEditFields();

	void Apollo7Padload();
	void Apollo8Padload();
	void Apollo9Padload();
	void Apollo10Padload();
	void Apollo11Padload(int LaunchDay);
	void Apollo12Padload();
	void Apollo13Padload();
	void Apollo14Padload();
	void Apollo15Padload();
	void Apollo16Padload();
	void Apollo17Padload();
public:
	CEdit LaunchMJDInput;
	afx_msg void OnBnClickedOk();
	afx_msg void OnCbnSelchangeCombo1();
	afx_msg void OnCbnSelchangeCombo2();
	afx_msg void OnCbnSelchangeCombo3();

	AGCPadloadGenerator agc;

	CComboBox Launchpad;
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
	CEdit LATSPLBox;
	CEdit LNGSPLBox;
	CEdit CSMMASSBox;
	CEdit LEMMASSBox;
	CEdit PACTOFFBox;
	CEdit YACTOFFBox;
	CEdit LADPADBox;
	CEdit LODPADBox;
	CEdit ALFAPADBox;
	CEdit P37RANGEBox;
	CEdit EMDOTBox;
	CEdit MinImp1Box;
	CEdit MinImp2Box;
	CEdit MinImp3Box;
	CEdit MinImp4Box;
	CStatic MinImp1Label;
	CStatic MinImp2Label;
	CStatic MinImp3Label;
	CStatic MinImp4Label;
	CStatic MinImp1Unit;
	CStatic MinImp2Unit;
	CStatic MinImp3Unit;
	CStatic MinImp4Unit;
	CEdit TRUNSFBox;
	CEdit SHAFTSFBox;
	CEdit WMIDPOSBox;
	CEdit WMIDVELBox;
	CEdit RVARMINBox;
	CEdit OutputBox;
	CEdit PIOSDataSetBox;
	CButton R2ModelBox;
	CEdit PBIASXBox;
	CEdit PIPASCFXBox;
	CEdit PBIASYBox;
	CEdit PIPASCFYBox;
	CEdit PBIASZBox;
	CEdit PIPASCFZBox;
	CEdit NBDXBox;
	CEdit NBDYBox;
	CEdit NBDZBox;
	CEdit ADIAXBox;
	CEdit ADIAYBox;
	CEdit ADIAZBox;
	CEdit ADSRAXBox;
	CEdit ADSRAYBox;
	CEdit ADSRAZBox;
};
