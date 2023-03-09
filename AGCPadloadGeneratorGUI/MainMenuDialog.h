#pragma once

class CAGCPadloadGeneratorGUIDlg;
class LGCPadloadGenerator;

class MainMenuDialog : public CDialogEx
{
	DECLARE_DYNAMIC(MainMenuDialog)

public:
	MainMenuDialog(CWnd* pParent = nullptr);   // Standardkonstruktor
	virtual ~MainMenuDialog();

// Dialogfelddaten
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_DIALOG1 };
#endif

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV-Unterstützung

	CAGCPadloadGeneratorGUIDlg *cmc_form;
	LGCPadloadGenerator *lgc_form;

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedCMC();
	afx_msg void OnBnClickedLGC();
};
