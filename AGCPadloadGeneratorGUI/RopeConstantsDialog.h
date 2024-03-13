#pragma once


// RopeConstantsDialog dialog

class RopeConstantsDialog : public CDialogEx
{
	DECLARE_DYNAMIC(RopeConstantsDialog)

public:
	RopeConstantsDialog(CWnd* pParent = nullptr);   // standard constructor
	virtual ~RopeConstantsDialog();

// Dialog Data
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_DIALOG5 };
#endif

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedGenerate();
	CEdit RopeConstantsYear;
	CEdit OutputBox;
};
